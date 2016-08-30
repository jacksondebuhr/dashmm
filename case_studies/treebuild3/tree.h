// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license. See the LICENSE file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================

#ifndef __TREE_H__
#define __TREE_H__


// C library
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cassert>

// C++ library
#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>

// HPX-5
#include "hpx/hpx.h"

// DASHMM
#include "dashmm/array.h"
#include "dashmm/dag.h"
#include "dashmm/index.h"
#include "dashmm/point.h"
#include "dashmm/domaingeometry.h"
#include "dashmm/reductionops.h"
#include "dashmm/shareddata.h"
using dashmm::int_sum_ident_op;
using dashmm::int_sum_op;
using dashmm::SharedData;
using dashmm::DomainGeometry;
using dashmm::Point;
using dashmm::Index;
using dashmm::Array;
using dashmm::ArrayRef;
using dashmm::ArrayData;
using dashmm::DAGInfo;
using dashmm::DAGNode;
using dashmm::DAGEdge;

// New things for this example
#include "rankwise.h"


// TODO - this will be removed eventually in favor of a foward declaration of
// the evaluator
template <typename Record>
class NodeRegistrar;

template <typename Source, typename Target, typename Record,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class Method,
          typename DistroPolicy>
class TreeRegistrar;

template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class Method,
          typename DistroPolicy>
class Registrar;


/// A Node of a tree.
///
/// This is a template over the record type. In DASHMM, this will be either
/// the Source type or the Target type.
template <typename Record>
class Node {
 public:
  using record_t = Record;
  using node_t = Node<Record>;
  using arrayref_t = ArrayRef<Record>;

  /// The default constructor allocates nothing, and sets all to zero.
  Node() : idx{}, parts{}, parent{nullptr}, first_{0} {
    for (int i = 0; i < 8; ++i) {
      child[i] = nullptr;
    }
    sema_ = HPX_NULL;
    complete_ = HPX_NULL;
  }

  /// Constuct with a known index. This will set the index of the node,
  /// and will allocate the completion detection LCO.
  Node(Index idx) : idx{idx}, parts{}, parent{nullptr}, first_{0} {
    for (int i = 0; i < 8; ++i) {
      child[i] = nullptr;
    }
    sema_ = HPX_NULL;
    complete_ = hpx_lco_and_new(8);
  }

  /// Construct with an index, a particle segment and a parent
  ///
  /// This will also create the completion detection LCO.
  Node(Index idx, arrayref_t parts, node_t *parent)
      : idx{idx}, parts{parts}, parent{parent}, first_{0} {
    for (int i = 0; i < 8; ++i) {
      child[i] = nullptr;
    }
    sema_ = HPX_NULL;
    complete_ = hpx_lco_and_new(8);
  }

  /// Any allocated LCOs are cleaned up explicitly.
  ///
  /// The completion LCOs are cleaned up as they finish their detection, and
  /// pass on any messages dependent on them. Very few nodes of the tree ever
  /// get a lock, so those are cleaned up when needed.
  ~Node() { }

  /// Returns the first record into which new records may be copied.
  ///
  /// NOTE: the node should be locked before using this to guarantee consistent
  /// access.
  size_t first() const {return first_;}

  /// Returns the LCO signaling completion of partitioning for this node
  hpx_addr_t complete() const {return complete_;}

  /// Gives the number of particles in the segment owned by this node
  size_t num_parts() const {return parts.n();}

  /// Create a completion detection and gate
  void add_completion() {
    assert(complete_ == HPX_NULL);
    complete_ = hpx_lco_and_new(8);
    assert(complete_ != HPX_NULL);
  }

  /// Create semaphore as needed
  void add_lock() {
    assert(sema_ == HPX_NULL);
    sema_ = hpx_lco_sema_new(1);
    assert(sema_ != HPX_NULL);
  }

  /// Delete semaphore
  void delete_lock() {
    assert(sema_ != HPX_NULL);
    hpx_lco_delete_sync(sema_);
    sema_ = HPX_NULL;
  }

  /// Lock the semaphore
  void lock() const {
    assert(sema_ != HPX_NULL);
    hpx_lco_sema_p(sema_);
  }

  /// Unlock the semaphore
  void unlock() const {
    assert(sema_ != HPX_NULL);
    hpx_lco_sema_v(sema_, HPX_NULL);
  }

  /// Increment the first record that is open
  ///
  /// See first() above. This will increment the location in which new records
  /// can be added to the rebalanced portion of the point data for this node.
  /// This will return if the given increment fills up the remaining spots in
  /// the segment of global array.
  ///
  /// NOTE: This should only be called once this node has been lock()-ed
  ///
  /// NOTE: This is only used for the nodes of the uniform grid.
  ///
  /// \param incr - the number of records that have just been added to the
  ///               segment.
  ///
  /// \returns - true if the segment is full; false otherwise
  bool increment_first(size_t incr) {
    first_ += incr;
    return first_ >= parts.n();
  }

  /// Partition the node
  ///
  /// In addition to sorting the points associated with this node, this will
  /// create the needed children and schedule the work of partitioning for
  /// those children.
  ///
  /// \param threshold - the partitioning threshold
  /// \param geo - the domain geometry
  // TODO update for same s and t
  void partition(int threshold, DomainGeometry *geo, int same_sandt) {
    size_t num_points = num_parts();
    assert(num_points >= 1);
    // TODO perhaps change the argument type to avoid this cast
    bool is_leaf = num_points <= (size_t)threshold;

    if (parent && parent->complete() != HPX_NULL) {
      hpx_call_when_with_continuation(complete_,
          parent->complete(), hpx_lco_set_action,
          complete_, hpx_lco_delete_action,
          nullptr, 0);
    }

    if (is_leaf) {
      // No children means this node is done with partitioning
      hpx_lco_and_set_num(complete_, 8, HPX_NULL);
    } else {
      // Compute center of the node
      double h = geo->size() / pow(2, idx.level());
      double center_x = geo->low().x() + (idx.x() + 0.5) * h;
      double center_y = geo->low().y() + (idx.y() + 0.5) * h;
      double center_z = geo->low().z() + (idx.z() + 0.5) * h;

      // Get the local data
      auto p_local = parts.pin();
      record_t *p = p_local.value();

      // Sort the particles among the children
      record_t *splits[9]{};
      splits[0] = p;
      splits[8] = &p[num_points];

      auto z_comp = [&center_z](record_t &a) {
        return a.position.z() < center_z;
      };
      splits[4] = std::partition(splits[0], splits[8], z_comp);

      auto y_comp = [&center_y](record_t &a) {
        return a.position.y() < center_y;
      };
      splits[2] = std::partition(splits[0], splits[4], y_comp);
      splits[6] = std::partition(splits[4], splits[8], y_comp);

      auto x_comp = [&center_x](record_t &a) {
        return a.position.x() < center_x;
      };
      splits[1] = std::partition(splits[0], splits[2], x_comp);
      splits[3] = std::partition(splits[2], splits[4], x_comp);
      splits[5] = std::partition(splits[4], splits[6], x_comp);
      splits[7] = std::partition(splits[6], splits[8], x_comp);

      // Perform some counting
      int stat[8]{};
      int offset[8]{};

      offset[0] = 0;
      stat[0] = splits[1] - splits[0];
      for (int i = 1; i < 8; ++i) {
        stat[i] = splits[i + 1] - splits[i];
        offset[i] = offset[i - 1] + stat[i - 1];
      }

      // Create child nodes
      for (int i = 0; i < 8; ++i) {
        if (stat[i]) {
          auto cparts = parts.slice(offset[i], stat[i]);
          node_t *cnd = new node_t{idx.child(i), cparts, this};
          child[i] = cnd;

          hpx_call(HPX_HERE, partition_node_, HPX_NULL,
                   &cnd, &geo, &threshold, &same_sandt);
        } else {
          hpx_lco_and_set(complete_, HPX_NULL);
        }
      }
    }
  }

  /// Return the size of the branch below this node
  ///
  /// This will return the total number of descendants of this node.
  ///
  /// \returns - the number of nodes in the branch below this node
  int n_descendants() const {
    // NOTE: the implementation is not recursive because HPX-5 can work with
    // small stack sizes.
    std::vector<const node_t *> V{this};
    int count = 0;
    while (!V.empty()) {
      count += V.size();
      std::vector<const node_t *> C{};
      for (size_t i = 0; i < V.size(); ++i) {
        for (int j = 0; j < 8; ++j) {
          if (V[i]->child[j]) {
            C.push_back(V[i]->child[j]);
          }
        }
      }
      V = std::move(C);
    }
    return count;
  }

  /// Return the number of immediate descendants of this node
  ///
  /// \returns - the number of non-null children of this node.
  int n_children() const {
    int retval{0};
    for (int i = 0; i < 8; ++i) {
      if (child[i] != nullptr) {
        ++retval;
      }
    }
    return retval;
  }

  /// Predicate for testing if this node is a leaf
  ///
  /// \returns - true if a leaf; false otherwise
  bool is_leaf() const {return n_children() == 0;}

  /// Compress the branch information into the provided buffers
  ///
  /// This will compress the branch information, which can then be sent
  /// to other ranks, where is can be extracted into the needed nodes.
  ///
  /// The structure of these is as two arrays, with one entry for each node
  /// except the root of the branch (which is the uniform grid node below
  /// which this branch is found, which is also this). One array gives the
  /// index of the parent of the current node, and the other gives which child
  /// of that parent is the current node.
  ///
  /// \param branch - which child of the parent indicated by tree is this node
  /// \param tree - which node is the parent of the given node
  void compress(int *branch, int *tree) const {
    // NOTE: the non-recursive implementation is to take into account HPX-5's
    // small default stack size.
    std::vector<const node_t *> V{this};
    std::vector<int> V_idx{-1};
    int curr = 0;
    while (!V.empty()) {
      std::vector<const node_t *> C{};
      std::vector<int> C_idx{};
      for (size_t i = 0; i < V.size(); ++i) {
        assert(V[i] != nullptr);
        for (int j = 0; j < 8; ++j) {
          if (V[i]->child[j]) {
            branch[curr] = j;
            tree[curr] = V_idx[i];
            C.push_back(V[i]->child[j]);
            C_idx.push_back(curr++);
          }
        }
      }
      V = std::move(C);
      V_idx = std::move(C_idx);
    }
  }

  /// Extract the branch information from the provided buffers
  ///
  /// This will extract the branch information received from other ranks,
  /// which can then be reified into the needed nodes.
  ///
  /// \param branch - buffers holding the tree structure
  /// \param tree -
  /// \param n_nodes - the number of nodes
  void extract(const int *branch, const int *tree, int n_nodes) {
    // Extract a compressed remote tree representation. As the tree is remote,
    // only {parent, child, idx} fields are needed.

    node_t *descendants = new node_t[n_nodes]{};

    // The compressed tree is created in depth first fashion. And there are two
    // choices here to fill in the parent, child, and idx fields of the
    // descendants.

    // Approach I: Setup parent and child, which is an embarassingly parallel
    // operation on the @p branch and @p tree. Afterwards, fan out along the
    // tree to fill in idx.

    // Approach II: Go over the input sequentially. For each node encountered,
    // by the depth first property, the index of its parent is already set.
    // So one can finish in one loop.

    // If on each rank, there are multiple subtrees being merged, approach II
    // might be sufficient. Approach I can be considered if finer granularity is
    // needed.

    // Approach II is implemented here.
    for (int i = 0; i < n_nodes; ++i) {
      int pos = tree[i];
      int which = branch[i];
      node_t *curr = &descendants[i];
      node_t *parent = (pos < 0 ? this : &descendants[pos]);

      curr->parent = parent;
      curr->idx = parent->idx.child(which);
      parent->child[which] = curr;
    }
  }

  /// Destroy the data in a branch
  ///
  /// This is only ever called on the local nodes. This will traverse the
  /// branch and delete all nodes below the given node.
  ///
  /// \param root - the root of the branch.
  static void destroy_branch(node_t *root) {
    // NOTE: this is not recursive because HPX-5 has small default stacks.
    std::vector<node_t *> V{root};
    while (!V.empty()) {
      std::vector<node_t *> C{};
      for (size_t i = 0; i < V.size(); ++i) {
        assert(V[i] != nullptr);
        for (int j = 0; j < 8; ++j) {
          if (V[i]->child[j]) {
            C.push_back(V[i]->child[j]);
          }
        }
        delete V[i];
        V[i] = nullptr;
      }
      V = std::move(C);
    }
  }


  Index idx;                      /// index of the node
  arrayref_t parts;               /// segment for this node
  node_t *parent;                 /// parent node
  node_t *child[8];               /// children of this node
  //DAGInfo dag;                    /// The DAG info for this node

  static hpx_action_t partition_node_;

 private:
  friend class NodeRegistrar<Record>;

  static int partition_node_handler(node_t *n, DomainGeometry *geo,
                                    int threshold, int same_sandt) {
    n->partition(threshold, geo, same_sandt);
    return HPX_SUCCESS;
  }

  size_t first_;            /// first record that is available
  hpx_addr_t sema_;         /// restrict concurrent modification
  hpx_addr_t complete_;     /// This is used to indicate that partitioning is
                            ///  complete
};

template <typename R>
hpx_action_t Node<R>::partition_node_ = HPX_ACTION_NULL;


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

// Tree needs a DualTree forward declaration
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class Method,
          typename DistroPolicy>
class DualTree;


/// The tree class represents one of the trees in a given evaluation
///
/// The Tree is a template over the typical DASHMM types, as well as a
/// Record type that indicates if this is a Source tree or a Target tree.
///
/// The nodes are arranged in a hybrid way in this tree. The top of the tree
/// (up to and including the finest uniform level) are allocated in an array.
/// The branches owned by this rank are allocated in the traditional fashion,
/// one node at a time. Branches from remotes are allocated in an array.
/// The two different schemes reflect that we sometimes know how many nodes we
/// shall need, and other times we do not.
template <typename Source, typename Target, typename Record,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class Method,
          typename DistroPolicy>
class Tree {
 public:
  using record_t = Record;
  using node_t = Node<Record>;
  using sourcenode_t = Node<Source>;
  using targetnode_t = Node<Target>;
  using arrayref_t = ArrayRef<Record>;
  using tree_t = Tree<Source, Target, Record, Expansion, Method, DistroPolicy>;
  using sourcetree_t = Tree<Source, Target, Source, Expansion, Method,
                            DistroPolicy>;
  using dualtree_t = DualTree<Source, Target, Expansion, Method, DistroPolicy>;

  /// Tree construction just default initializes the object
  Tree() : root_{nullptr}, unif_grid_{nullptr}, unif_done_{HPX_NULL},
           sorted_{} { }

  arrayref_t sorted() const {return sorted_;}

  /// Setup some basic information during initial tree construction
  ///
  /// This action is the target of a broadcast. The basic information about
  /// the tree is set up in this action. The most important of which is the
  /// allocation of the top portions of the tree containing the root, and
  /// extending to the finest uniform level of partitioning.
  static int setup_basics_handler(tree_t *tree, int unif_level) {
    tree->unif_done_ = hpx_lco_and_new(1);
    assert(tree->unif_done_ != HPX_NULL);

    // Setup unif_grid
    int n_top_nodes{1};
    for (int i = 1; i < unif_level; ++i) {
      n_top_nodes += pow(8, i);
    }
    int dim3 = pow(8, unif_level);
    tree->root_ = new node_t[n_top_nodes + dim3]{};
    tree->unif_grid_ = &tree->root_[n_top_nodes];

    tree->root_[0].idx = Index{0, 0, 0, 0};
    int startingnode{0};
    int stoppingnode{1};
    for (int level = 0; level < unif_level ; ++level) {
      for (int nd = startingnode; nd < stoppingnode; ++nd) {
        node_t *snode = &tree->root_[nd];

        if (level != unif_level - 1) {
          // At the lower levels, we do not require the ordering
          int firstchild = nd * 8 + 1;
          for (int i = 0; i < 8; ++i) {
            node_t *scnode = &tree->root_[firstchild + i];
            snode->child[i] = scnode;
            scnode->idx = snode->idx.child(i);
            scnode->parent = snode;
          }
        } else {
          // The final level does need the ordering
          for (int i = 0; i < 8; ++i) {
            Index cindex = snode->idx.child(i);
            uint64_t morton = morton_key(cindex.x(), cindex.y(), cindex.z());
            snode->child[i] = &tree->unif_grid_[morton];
            tree->unif_grid_[morton].idx = cindex;
            tree->unif_grid_[morton].parent = snode;
            tree->unif_grid_[morton].add_lock();
            tree->unif_grid_[morton].add_completion();
          }
        }
      }

      startingnode += pow(8, level);
      stoppingnode += pow(8, level + 1);
    }

    return HPX_SUCCESS;
  }

  /// Destroy allocated data for this tree.
  ///
  /// This will delete the local branches of the tree as well as destroying
  /// any locks allocated for the uniform grid, and will free the remote
  /// branches.
  static int delete_tree_handler(tree_t *tree, int ndim, int first, int last) {
    // NOTE: The difference here is that there are two different allocation
    // schemes for the nodes.

    for (int i = 0; i < first; ++i) {
      node_t *curr = &tree->unif_grid_[i];
      curr->delete_lock();

      for (int j = 0; j < 8; ++j) {
        node_t *child = curr->child[j];
        if (child) {
          delete [] child;
          break;
        }
      }
    }

    for (int i = first; i <= last; ++i) {
      node_t *curr = &tree->unif_grid_[i];
      curr->delete_lock();

      for (int j = 0; j < 8; ++j) {
        if (curr->child[j]) {
          node_t::destroy_branch(curr->child[j]);
        }
      }

      hpx_lco_delete_sync(curr->complete());
    }

    for (int i = last + 1; i < ndim; ++i) {
      node_t *curr = &tree->unif_grid_[i];
      curr->delete_lock();

      for (int j = 0; j < 8; ++j) {
        node_t *child = curr->child[j];
        if (child) {
          delete [] child;
          break;
        }
      }
    }

    delete [] tree->root_;

    hpx_lco_delete_sync(tree->unif_done_);

    return HPX_SUCCESS;
  }

  /// Assign points to the uniform grid.
  ///
  /// This will assign the points to the uniform grid. This gives the points
  /// the id (in the Morton Key sense) of the box to which they are assigned,
  /// and it will count the numbers in each box. This is a key first step to
  /// computing the distribution of the sources and targets during tree
  /// construction.
  static int assign_points_to_unif_grid(const record_t *P, int npts,
                                         const DomainGeometry *geo,
                                         int unif_level, int *gid,
                                         int *count) {
    Point corner = geo->low();
    double scale = 1.0 / geo->size();

    // TODO: This is serial processing; is there some way to parallelize this?
    //   This would perhaps be worth timing.
    int dim = pow(2, unif_level);
    for (int i = 0; i < npts; ++i) {
      const record_t *p = &P[i];
      int xid = std::min(dim - 1,
                         (int)(dim * (p->position.x() - corner.x()) * scale));
      int yid = std::min(dim - 1,
                         (int)(dim * (p->position.y() - corner.y()) * scale));
      int zid = std::min(dim - 1,
                         (int)(dim * (p->position.z() - corner.z()) * scale));
      gid[i] = morton_key(xid, yid, zid);
      count[gid[i]]++;
    }

    return HPX_SUCCESS;
  }

  /// Reorder the particles according to their place in the uniform grid
  ///
  /// This will rearrange the particles into their bin order. This is a stable
  /// reordering.
  ///
  /// This routine is adapted from a routine in the publicly available code
  /// GADGET-2 (http://wwwmpa.mpa-garching.mpg.de/gadget/).
  //
  // TODO: Note that this too is a serial operation. Could we do some on-rank
  // parallelism here?
  static int *group_points_on_unif_grid(record_t *p_in, int npts,
                                        int dim3, int *gid_of_points,
                                        const int *count, int **retval) {
    int *offset = new int[dim3]();

    offset[0] = 0;
    for (int i = 1; i < dim3; ++i) {
      offset[i] = offset[i - 1] + count[i - 1];
    }

    // Set gid to the final location of the particle; NOTE: this will modify
    // the offset. This is corrected below
    for (int i = 0; i < npts; ++i) {
      int gid = gid_of_points[i];
      gid_of_points[i] = offset[gid];
      offset[gid]++;
    }

    // Declare some loop variables
    record_t source, save;
    int source_sort, save_sort;
    int isource, isave, idest;

    // Do an O(N) rearrangement
    for (int i = 0; i < npts; ++i) {
      if (gid_of_points[i] != i) {
        source = p_in[i];
        source_sort = gid_of_points[i];
        isource = gid_of_points[i];
        idest = gid_of_points[i];

        do {
          save = p_in[idest];
          save_sort = gid_of_points[idest];
          isave = gid_of_points[idest];

          p_in[idest] = source;
          gid_of_points[idest] = source_sort;

          if (idest == i) break;

          source = save;
          source_sort = save_sort;
          isource = isave;

          idest = isource;
        } while (1);
      }
    }

    // Correct offset
    for (int i = 0; i < dim3; ++i) {
      offset[i] -= count[i];
    }

    *retval = offset;

    return HPX_SUCCESS;
  }

  /// Receive partitioned tree nodes from a remote locality
  ///
  /// This action receives the compressed representation of a remote branch
  /// of the tree, and then expands it into node objects, connectint it to
  /// the correct portion of the tree at this locality. In general, there will
  /// be one such message per uniform grid node not owned by this locality, per
  /// tree.
  static int recv_node_handler(char *message_buffer, size_t UNUSED) {
    hpx_addr_t *rwdata = reinterpret_cast<hpx_addr_t *>(message_buffer);
    int *compressed_tree = reinterpret_cast<int *>(
                                message_buffer + sizeof(hpx_addr_t));
    RankWise<dualtree_t> global_tree{*rwdata};
    auto local_tree = global_tree.here();
    int id = compressed_tree[0];
    int n_nodes = compressed_tree[1];
    int type = compressed_tree[2];

    // NOTE: This is a miserable solution to this problem. See if there is
    // something better for this.
    if (type == 0) {
      Node<Source> *grid = local_tree->unif_grid_source();

      if (n_nodes) {
        const int *branch = &compressed_tree[3];
        const int *tree = &compressed_tree[3 + n_nodes];
        grid[id].extract(branch, tree, n_nodes);
      }

      hpx_lco_and_set_num(grid[id].complete(), 8, HPX_NULL);
    } else {
      Node<Target> *grid = local_tree->unif_grid_target();

      if (n_nodes) {
        const int *branch = &compressed_tree[3];
        const int *tree = &compressed_tree[3 + n_nodes];
        grid[id].extract(branch, tree, n_nodes);
      }

      hpx_lco_and_set_num(grid[id].complete(), 8, HPX_NULL);
    }

    return HPX_SUCCESS;
  }

  /// This action sends compressed node representation to remote localities.
  ///
  /// This action is called once the individual grids are done.
  static int send_node_handler(node_t *n, record_t *sorted,
                               int id, hpx_addr_t rwaddr, int type) {
    node_t *curr = &n[id];

    // Exclude curr as it is already allocated on remote localities
    int n_nodes = curr->n_descendants() - 1;
    size_t msgsize = sizeof(int) * (3 + n_nodes * 2) + sizeof(hpx_addr_t);
    char *message_buffer = new char[msgsize];
    hpx_addr_t *rwdata = reinterpret_cast<hpx_addr_t *>(message_buffer);
    int *compressed_tree = reinterpret_cast<int *>(
                                message_buffer + sizeof(hpx_addr_t));

    *rwdata = rwaddr;
    compressed_tree[0] = id; // where to merge
    compressed_tree[1] = n_nodes; // # of nodes
    compressed_tree[2] = type; //type
    if (n_nodes) {
      int *branch = &compressed_tree[3];
      int *tree = &compressed_tree[3 + n_nodes];
      curr->compress(branch, tree);
    }

    int rank = hpx_get_my_rank();
    int num_ranks = hpx_get_num_ranks();
    hpx_addr_t done = hpx_lco_and_new(num_ranks - 1);

    for (int r = 0; r < num_ranks; ++r) {
      if (r != rank) {
        hpx_parcel_t *p = hpx_parcel_acquire(message_buffer, msgsize);
        hpx_parcel_set_target(p, HPX_THERE(r));
        hpx_parcel_set_action(p, recv_node_);
        hpx_parcel_send(p, done);
      }
    }

    // Clear local memory once all parcels are sent
    hpx_lco_wait(done);

    delete [] message_buffer;
    return HPX_SUCCESS;
  }

  /// Initialize point exchange between ranks
  ///
  /// This is largely just a setup procedure to compute a few things that are
  /// needed for point exchanges. The most important of which is the global
  /// offset in this rank's data where each incoming batch of points will be
  /// stored in the sorted array. Also, this allocates the segment of memory
  /// that will store the sorted data.
  ///
  /// \param tree - the tree object
  /// \param first - the first node of the uniform grid owned by this rank
  /// \param last - the last node of the uniform grid owned by this rank
  /// \param global_count - the global counts in each node of the uniform grid
  /// \param local_count - the local counts in each node of the uniform grid
  /// \param local_offset - where in the local data are the points for each node
  ///                       of the uniform grid
  /// \param geo - the domain geometry for the tree
  /// \param threshold - the partitioning threshold for the tree
  /// \param temp - the local point data
  /// \param n - the uniform grid nodes
  ///
  /// \returns - the global offset into this rank's data for each uniform grid
  ///            node's points.
  static void init_point_exchange(tree_t *tree, int first, int last,
                                  const int *global_count,
                                  const int *local_count,
                                  const int *local_offset,
                                  const DomainGeometry *geo,
                                  int threshold,
                                  const record_t *temp,
                                  node_t *n) {
    int range = last - first + 1;
    arrayref_t sorted_ref{};

    // Compute global_offset
    int *global_offset = new int[range]();
    size_t num_points = global_count[first];
    for (int i = first + 1; i <= last; ++i) {
      num_points += global_count[i];
      global_offset[i - first] = global_offset[i - first - 1] +
                                 global_count[i - 1];
    }

    if (num_points > 0) {
      hpx_addr_t sorted_gas = hpx_gas_alloc_local(
          1, sizeof(record_t) * num_points, 0);
      assert(sorted_gas != HPX_NULL);
      sorted_ref = arrayref_t{sorted_gas, num_points, num_points};

      for (int i = first; i <= last; ++i) {
        node_t *curr = &n[i];
        curr->parts = sorted_ref.slice(global_offset[i - first],
                                       global_count[i]);

        if (local_count[i]) {
          // Copy local points before merging remote points
          curr->lock();
          auto sorted = curr->parts.pin();
          memcpy(sorted.value() + curr->first(),
                 &temp[local_offset[i]],
                 sizeof(record_t) * local_count[i]);

          if (curr->increment_first(local_count[i])) {
            // This grid does not expect remote points.
            // Spawn adaptive partitioning
            int ssat = 0;
            hpx_call(HPX_HERE, node_t::partition_node_, HPX_NULL,
                     &curr, &geo, &threshold, &ssat);
          }
          curr->unlock();
        }
      }
    }

    // TODO this will have to be different - the target tree might need to use
    // the one from the source tree
    tree->sorted_ = sorted_ref;

    hpx_lco_and_set(tree->unif_done_, HPX_NULL);

    delete [] global_offset;
  }

  static void init_point_exchange_same_s_and_t(
      tree_t *tree, int first, int last, const int *global_count,
      const int *local_count, const int *local_offset,
      const DomainGeometry *geo, int threshold, const record_t *temp,
      node_t *n, sourcenode_t *snodes, sourcetree_t *source_tree) {
    int range = last - first + 1;
    arrayref_t sorted_ref{};

    // Compute global_offset
    int *global_offset = new int[range]();
    size_t num_points = global_count[first];
    for (int i = first + 1; i <= last; ++i) {
      num_points += global_count[i];
      global_offset[i - first] = global_offset[i - first - 1] +
                                 global_count[i - 1];
    }

    if (num_points > 0) {
      for (int i = first; i <= last; ++i) {
        // When S == T, we need to set the parts on the target to be the
        // parts on the equivalent source node
        node_t *curr = &n[i];
        sourcenode_t *curr_source = &snodes[i];
        ArrayRef<Source> sparts = curr_source->parts;
        curr->parts = arrayref_t{sparts.data(), sparts.n(), sparts.n_tot()};

        if (local_count[i]) {
          curr->lock();
          // still need to increment first, and still need to partition
          // if things are ready. The only catch is that we have to wait for
          // the equivalent source node to be ready before we start.
          if (curr->increment_first(local_count[i])) {
            // This grid does not expect remote points.
            // Spawn adaptive partitioning
            int ssat = 1;
            hpx_call_when(curr_source->complete(), HPX_HERE,
                          node_t::partition_node_, HPX_NULL,
                          &curr, &geo, &threshold, &ssat);
          }
          curr->unlock();
        }
      }
    }

    // Again, the target tree reuses the data from the source tree
    ArrayRef<Source> ssort = source_tree->sorted();
    tree->sorted_ = arrayref_t{ssort.data(), ssort.n(), ssort.n_tot()};

    hpx_lco_and_set(tree->unif_done_, HPX_NULL);

    delete [] global_offset;
  }

  /// Merge incoming points into the local array
  ///
  /// This action merges particular points with the sorted list. Also, if this
  /// is the last set of points that are merged, this will go ahead and start
  /// the adaptive partitioning of that part of the local tree.
  ///
  /// \param temp - the local data
  /// \param n - the uniform grid node
  /// \param n_arrived - the number arriving in this message
  /// \param rwgas - the address of the rankwise dual tree
  static int merge_points_handler(record_t *temp, node_t *n, int n_arrived,
                                  hpx_addr_t rwgas) {
    RankWise<dualtree_t> global_tree{rwgas};
    auto local_tree = global_tree.here();

    // Note: all the pointers are local to the calling rank.
    n->lock();
    size_t first = n->first();
    auto localp = n->parts.pin();
    record_t *p = localp.value();
    memcpy(p + first, temp, sizeof(record_t) * n_arrived);

    if (n->increment_first(n_arrived)) {
      const DomainGeometry *geoarg = local_tree->domain();
      int thresh = local_tree->refinement_limit();
      int ssat = 0;
      hpx_call(HPX_HERE, node_t::partition_node_, HPX_NULL,
               &n, &geoarg, &thresh, &ssat);
    }
    n->unlock();

    return HPX_SUCCESS;
  }

  /// Merge incoming points into the local array
  ///
  /// This action performs the S == T version of point merging. In this version,
  /// no points are actually merged, as that is handles in the Source version of
  /// this routine. Instead, this merely updates the first counter, and
  /// calls to partitioning when appropriate. Note, this is now a call when
  /// waiting on the completion of the partitioning in the source tree.
  ///
  /// \param temp - the local data
  /// \param n - the uniform grid node
  /// \param n_arrived - the number arriving in this message
  /// \param rwgas - the address of the rankwise dual tree
  // TODO very little goes on in this. Perhaps hoist this into the calling
  // context.
  static int merge_points_same_s_and_t_handler(targetnode_t *target_node,
      int n_arrived, sourcenode_t *source_node, hpx_addr_t rwgas) {
    RankWise<dualtree_t> global_tree{rwgas};
    auto local_tree = global_tree.here();

    target_node->lock();
    if (target_node->increment_first(n_arrived)) {
      const DomainGeometry *geoarg = local_tree->domain();
      int thresh = local_tree->refinement_limit();
      int ssat = 1;
      hpx_call_when(source_node->complete(),
                    HPX_HERE, node_t::partition_node_, HPX_NULL,
                    &target_node, &geoarg, &thresh, &ssat);
    }
    target_node->unlock();

    return HPX_SUCCESS;
  }

  /// Split the bits of an integer to be used in a Morton Key
  static uint64_t split(unsigned k) {
    uint64_t split = k & 0x1fffff;
    split = (split | split << 32) & 0x1f00000000ffff;
    split = (split | split << 16) & 0x1f0000ff0000ff;
    split = (split | split << 8)  & 0x100f00f00f00f00f;
    split = (split | split << 4)  & 0x10c30c30c30c30c3;
    split = (split | split << 2)  & 0x1249249249249249;
    return split;
  }

  /// Compute the Morton key for a gives set of indices
  static uint64_t morton_key(unsigned x, unsigned y, unsigned z) {
    uint64_t key = 0;
    key |= split(x) | split(y) << 1 | split(z) << 2;
    return key;
  }


private:
  friend class DualTree<Source, Target, Expansion, Method, DistroPolicy>;
  friend class TreeRegistrar<Source, Target, Record, Expansion, Method,
                             DistroPolicy>;

  node_t *root_;            /// Root of the tree
  node_t *unif_grid_;       /// The uniform grid
  hpx_addr_t unif_done_;    /// An LCO indicating that the uniform partition is
                            /// complete
  arrayref_t sorted_;       /// A reference to the sorted point data owned by
                            /// this tree.

  static hpx_action_t setup_basics_;
  static hpx_action_t delete_tree_;
  static hpx_action_t recv_node_;
  static hpx_action_t send_node_;
  static hpx_action_t assign_points_;
  static hpx_action_t group_points_;
  static hpx_action_t merge_points_;
  static hpx_action_t merge_points_same_s_and_t_;
};

template <typename S, typename T, typename R,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, R, E, M, D>::setup_basics_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, R, E, M, D>::delete_tree_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, R, E, M, D>::recv_node_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, R, E, M, D>::send_node_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, R, E, M, D>::assign_points_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, R, E, M, D>::group_points_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, R, E, M, D>::merge_points_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, R, E, M, D>::merge_points_same_s_and_t_
    = HPX_ACTION_NULL;


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


/// The DualTree organizes the source and target tree and handles common work
///
/// The DualTree manages all work that instersects between the two trees.
/// This object stores the pointers to the source and target trees.
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class Method,
          typename DistroPolicy>
class DualTree {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Method<Source, Target, Expansion, DistroPolicy>;
  using sourcenode_t = Node<Source>;
  using targetnode_t = Node<Target>;
  using sourceref_t = ArrayRef<Source>;
  using sourcearraydata_t = ArrayData<Source>;
  using targetref_t = ArrayRef<Target>;
  using targetarraydata_t = ArrayData<Target>;
  using sourcetree_t = Tree<Source, Target, Source, Expansion, Method,
                            DistroPolicy>;
  using targettree_t = Tree<Source, Target, Target, Expansion, Method,
                            DistroPolicy>;
  using dualtree_t = DualTree<Source, Target, Expansion, Method, DistroPolicy>;

  /// Construction is always default
  DualTree()
    : domain_{}, refinement_limit_{1}, unif_level_{1}, dim3_{8},
      unif_count_{HPX_NULL}, unif_count_value_{nullptr},
      distribute_{nullptr}, method_{}, source_tree_{nullptr},
      target_tree_{nullptr} { }

  /// We delete the copy constructor and copy assignement operator.
  DualTree(const dualtree_t &other) = delete;
  dualtree_t &operator=(const dualtree_t &other) = delete;

  // TODO define move construction and assignment?
  // Is this needed?

  /// Return the number of source points in the uniform grid
  ///
  /// This routine returns the address of the given information,
  ///
  /// \parma i - the uniform grid node in question
  int *unif_count_src(size_t i = 0) const {return &unif_count_value_[i];}

  /// Return the number of target points in the uniform grid
  ///
  /// This routine returns the address of the given information,
  ///
  /// \parma i - the uniform grid node in question
  int *unif_count_tar(size_t i = 0) const {
    return &unif_count_value_[i + dim3_];
  }

  /// Return the DomainGeometry for this DualTree.
  const DomainGeometry *domain() const {return &domain_;}

  /// Return the refinement limit used to build the tree.
  int refinement_limit() const {return refinement_limit_;}

  /// Return the uniform grid for the source tree.
  sourcenode_t *unif_grid_source() {return source_tree_->unif_grid_;}

  /// Return the uniform grid for the target tree.
  targetnode_t *unif_grid_target() {return target_tree_->unif_grid_;}

  /// Return the number of post-sorting sources.
  size_t sorted_src_count() const {return source_tree_->sorted_.n_tot();}

  /// Return the number of post-sorting targets.
  size_t sorted_tar_count() const {return target_tree_->sorted_.n_tot();}

  /// Return the sorted sources.
  sourcearraydata_t sorted_src() const {return source_tree_->sorted_.pin();}

  /// Return the sorted targets.
  targetarraydata_t sorted_tar() const {return target_tree_->sorted_.pin();}

  /// Return the method this object will use for DAG operations.
  const method_t &method() const {return method_;}

  /// Set the method this object will use for DAG operations.
  void set_method(const method_t &method) {method_ = method;}

  /// Destroy any allocated memory associated with this DualTree.
  void clear_data() {
    int rank = hpx_get_my_rank();

    int b = first(rank);
    int e = last(rank);

    // Tell each tree to delete itself
    hpx_addr_t clear_done = hpx_lco_and_new(2);
    assert(clear_done != HPX_NULL);
    hpx_call(HPX_HERE, sourcetree_t::delete_tree_, clear_done,
             &source_tree_, &dim3_, &b, &e);
    hpx_call(HPX_HERE, targettree_t::delete_tree_, clear_done,
             &target_tree_, &dim3_, &b, &e);
    hpx_lco_wait(clear_done);
    hpx_lco_delete_sync(clear_done);

    // the unif counts
    delete [] unif_count_value_;

    // then free the trees themselves
    delete source_tree_;
    delete target_tree_;

    delete [] distribute_;
  }

  /// Return the first uniform grid node owned by the given rank.
  int first(int rank) const {return rank == 0 ? 0 : distribute_[rank - 1] + 1;}

  /// Return the last uniform grid node owned by the given rank.
  int last(int rank) const {return distribute_[rank];}


  // More DASHMM stuff
  // DAG create_DAG(bool same_sandt)
  // void collect_DAG_nodes(DAG &dag)
  // void create_expansions_from_DAG(int n_digits)
  // hpx_addr_t setup_termination_detection(DAG &dag)
  // void setup_edge_lists(DAG &dag)
  // void start_DAG_evaluation()
  // void destroy_DAG_LCOs(DAG &dag)


  /// Create the basic data for a distributed tree for use with DASHMM
  ///
  /// This will compute the domain from the given source and target points and
  /// will create the basic data for a distributed tree with that domain. The
  /// return object is a RankWise object storing each copy of the local tree.
  ///
  /// This call is synchronous, and should be called from inside an HPX thread.
  /// Further, this is to be called in a diffusive style; only a single thread
  /// should call this function.
  ///
  /// \param threshold - the partitioning threshold for the tree
  /// \param sources - the source data
  /// \param targets - the target data
  ///
  /// \returns - the RankWise object containing the dual tree
  static RankWise<dualtree_t> create(int threshold, Array<Source> sources,
                                     Array<Target> targets) {
    bool same_sandt{false};
    if (sources.data() == targets.data()) {
      same_sandt = true;
    }
    hpx_addr_t domain_geometry = compute_domain_geometry(sources, targets,
                                                         same_sandt);
    RankWise<dualtree_t> retval = setup_basic_data(threshold, domain_geometry,
                                                   same_sandt);
    hpx_lco_delete_sync(domain_geometry);
    return retval;
  }

  /// Partition the tree
  ///
  /// This will do the bulk of the work for partitioning and creating the
  /// distributed tree. After the call to this routine, the source and target
  /// data will be redistributed to make evaluation easier for DASHMM. Behind
  /// the scenes, the local segments of the arrays inside sources and targets
  /// will have been replaced.
  ///
  /// This routine should be called from inside an HPX thread. Further, this is
  /// to be called in a diffusive style; this should be called from only a
  /// single thread. This routine will handle involving all ranks in the
  /// system.
  ///
  /// \param global_tree - an object previously initialized with create()
  /// \param sources - the source data
  /// \param targets - the target data
  ///
  /// \returns - an LCO indication completion of the partitioning.
  static hpx_addr_t partition(RankWise<dualtree_t> global_tree,
                              Array<Source> sources,
                              Array<Target> targets) {
    hpx_addr_t retval = hpx_lco_future_new(0);
    assert(retval != HPX_NULL);

    hpx_addr_t tree_gas = global_tree.data();
    hpx_addr_t source_gas = sources.data();
    hpx_addr_t target_gas = targets.data();
    hpx_bcast_lsync(create_dual_tree_, retval,
                    &tree_gas, &source_gas, &target_gas);

    return retval;
  }

  /// Destroy a distributed tree.
  ///
  /// This cleans up all allocated resources used by the DualTree.
  ///
  /// This should be called from an HPX thread, in a diffusive style.
  ///
  /// \param global_tree - the distributed tree
  static void destroy(RankWise<dualtree_t> global_tree) {
    hpx_addr_t rwtree = global_tree.data();
    hpx_bcast_rsync(finalize_partition_, &rwtree);

    auto tree = global_tree.here();
    hpx_lco_delete_sync(tree->unif_count_);
  }


 private:
  friend class Registrar<Source, Target, Expansion, Method, DistroPolicy>;

  /// Action to set the domain geometry given the sources and targets
  ///
  /// This action is the target of a broadcast, and computes the domain for
  /// the local sources and targets. The result is then given to a reduction
  /// LCO which reduces each rank's portion.
  ///
  /// \param sources_gas - the global address of the source data
  /// \param targets_gas - the global address of the target data
  /// \param domain_geometry - a reduction LCO to which the local domain is sent
  static int set_domain_geometry_handler(hpx_addr_t sources_gas,
                                         hpx_addr_t targets_gas,
                                         hpx_addr_t domain_geometry,
                                         int same_sandt) {
    Array<source_t> sources{sources_gas};
    sourceref_t src_ref = sources.ref();
    sourcearraydata_t src_data = src_ref.pin();
    Source *s = src_data.value();

    double var[6] = {1e50, -1e50, 1e50, -1e50, 1e50, -1e50};

    // TODO: add more parallelism
    for (size_t i = 0; i < src_ref.n(); ++i) {
      var[0] = fmin(var[0], s[i].position.x());
      var[1] = fmax(var[1], s[i].position.x());
      var[2] = fmin(var[2], s[i].position.y());
      var[3] = fmax(var[3], s[i].position.y());
      var[4] = fmin(var[4], s[i].position.z());
      var[5] = fmax(var[5], s[i].position.z());
    }

    // Only do the targets if they are different from the sources
    if (!same_sandt) {
      Array<target_t> targets{targets_gas};
      targetref_t trg_ref = targets.ref();
      targetarraydata_t trg_data = trg_ref.pin();
      Target *t = trg_data.value();

      for (size_t i = 0; i < trg_ref.n(); ++i) {
        var[0] = fmin(var[0], t[i].position.x());
        var[1] = fmax(var[1], t[i].position.x());
        var[2] = fmin(var[2], t[i].position.y());
        var[3] = fmax(var[3], t[i].position.y());
        var[4] = fmin(var[4], t[i].position.z());
        var[5] = fmax(var[5], t[i].position.z());
      }
    }

    hpx_lco_set_lsync(domain_geometry, sizeof(double) * 6, var, HPX_NULL);

    return HPX_SUCCESS;
  }

  /// Operation implementing the identity for the domain reduction
  static void domain_geometry_init_handler(double *values,
                                           const size_t UNUSED) {
    values[0] = 1e50; // xmin
    values[1] = -1e50; // xmax
    values[2] = 1e50; // ymin
    values[3] = -1e50; // ymax
    values[4] = 1e50; // zmin
    values[5] = -1e50; // zmax
  }

  /// Operation implementing the reduction for the domain reduction
  static void domain_geometry_op_handler(double *lhs, double *rhs,
                                         size_t UNUSED) {
    lhs[0] = fmin(lhs[0], rhs[0]);
    lhs[1] = fmax(lhs[1], rhs[1]);
    lhs[2] = fmin(lhs[2], rhs[2]);
    lhs[3] = fmax(lhs[3], rhs[3]);
    lhs[4] = fmin(lhs[4], rhs[4]);
    lhs[5] = fmax(lhs[5], rhs[5]);
  }

  /// Compute the bounding box for the given source and target points.
  ///
  /// This will return an LCO which will contain six doubles, (xmin, xmax,
  /// ymin, ymax, zmin, zmax) This is an asynchronous call, it will return
  /// as soon as the work is scheduled.
  ///
  /// \param sources - the source array
  /// \param targets - the target array
  ///
  /// \returns - address of an LCO containing the reduced domain
  static hpx_addr_t compute_domain_geometry(Array<Source> sources,
                                            Array<Target> targets,
                                            bool same_sandt) {
    // Create a reduction LCO
    hpx_addr_t domain_geometry =
      hpx_lco_reduce_new(hpx_get_num_ranks(), sizeof(double) * 6,
                         domain_geometry_init_,
                         domain_geometry_op_);

    // Launch the reduction actions
    hpx_addr_t sglob = sources.data();
    hpx_addr_t tglob = targets.data();
    int ssat = (same_sandt ? 1 : 0);
    hpx_bcast_lsync(set_domain_geometry_, HPX_NULL,
                    &sglob, &tglob, &domain_geometry, &ssat);

    return domain_geometry;
  }

  /// Action to perform initializtion of basic data for the local tree
  ///
  /// This is the target of a broadcast, and it sets various data about the
  /// tree.
  ///
  /// \param rwdata - the global address of the global tree
  /// \param count - an LCO in which the uniform grid counting is reduced
  /// \param limit - the partitioning threshold for the tree
  /// \param domain_geometry - the LCO in which the domain is reduced
  static int init_partition_handler(hpx_addr_t rwdata, hpx_addr_t count,
                                    int limit, hpx_addr_t domain_geometry,
                                    int same_sandt) {
    RankWise<dualtree_t> global_tree{rwdata};
    auto tree = global_tree.here();

    int num_ranks = hpx_get_num_ranks();
    tree->unif_level_ = ceil(log(num_ranks) / log(8)) + 1;
    tree->dim3_ = pow(8, tree->unif_level_);
    tree->unif_count_ = count;
    tree->refinement_limit_ = limit;
    tree->source_tree_ =
      new Tree<Source, Target, Source, Expansion, Method, DistroPolicy>{};
    tree->target_tree_ =
      new Tree<Source, Target, Target, Expansion, Method, DistroPolicy>{};
    tree->same_sandt_ = same_sandt;

    // Call out to tree setup stuff
    hpx_addr_t setup_done = hpx_lco_and_new(2);
    assert(setup_done != HPX_NULL);

    hpx_call(HPX_HERE, sourcetree_t::setup_basics_, setup_done,
             &tree->source_tree_, &tree->unif_level_);
    hpx_call(HPX_HERE, targettree_t::setup_basics_, setup_done,
             &tree->target_tree_, &tree->unif_level_);

    // We here allocate space for the result of the counting
    tree->unif_count_value_ = new int[tree->dim3_ * 2]();

    // Setup domain_
    double var[6];
    hpx_lco_get(domain_geometry, sizeof(double) * 6, &var);
    double length = fmax(var[1] - var[0],
                         fmax(var[3] - var[2], var[5] - var[4]));
    DomainGeometry geo{Point{(var[1] + var[0] - length) / 2,
                             (var[3] + var[2] - length) / 2,
                             (var[5] + var[4] - length) / 2}, length};
    tree->domain_ = geo;

    // Wait for setup to be done
    hpx_lco_wait(setup_done);
    hpx_lco_delete_sync(setup_done);

    return HPX_SUCCESS;
  }

  /// Allocate and setup a Dual Tree
  ///
  /// This will both allocate and setup a dual tree.
  ///
  /// \param threshold - the partitioning threshold
  /// \param domain_geometry - an LCO into which the domain is reduced
  static RankWise<dualtree_t> setup_basic_data(int threshold,
                                               hpx_addr_t domain_geometry,
                                               bool same_sandt) {
    RankWise<dualtree_t> retval{};
    retval.allocate();
    if (!retval.valid()) {
      // We return the invalid value to indicate the error
      return retval;
    }

    // Now the single things are created.
    int num_ranks = hpx_get_num_ranks();
    int level = ceil(log(num_ranks) / log(8)) + 1;
    int dim3 = pow(8, level);
    hpx_addr_t ucount = hpx_lco_reduce_new(num_ranks, sizeof(int) * (dim3 * 2),
                                           int_sum_ident_op,
                                           int_sum_op);
    hpx_addr_t rwdata = retval.data();
    int ssat = (same_sandt ? 1 : 0);
    hpx_bcast_rsync(init_partition_, &rwdata, &ucount, &threshold,
                    &domain_geometry, &ssat);

    return retval;
  }

  // TODO: Note that this would be a target for a second member of the
  // Distribution Policy
  // TODO: can this be simplified? In particular, use the same target for each
  // rather than a changing target.
  //
  /// Partition the available points among the localities
  ///
  /// Given the uniform counts, this will provide a good guess at a
  /// distribution of those points. This operates basically through the Morton
  /// ordering of the uniform grid nodes, and aims to have segments of the
  /// space-filling curve which have similar total counts.
  ///
  /// \param num_ranks - the number of localities to divide between
  /// \param global - the global counts
  /// \param len - the number of unform grid nodes
  static int *distribute_points(int num_ranks, const int *global, int len) {
    int *ret = new int[num_ranks]();

    const int *s = global; // Source counts
    const int *t = &global[len]; // Target counts

    int total = 0;
    for (int i = 0; i < len; ++i) {
      total += s[i] + t[i];
    }

    int rank = 0;
    int iterator = 0;

    while (rank < num_ranks && total > 0) {
      if (rank == num_ranks - 1) {
        // Take the remaining grids
        ret[rank++] = len - 1;
      } else {
        int avg = total / (num_ranks - rank);
        int sum = 0;
        int sum1 = 0;

        for (int i = iterator; i < len; ++i) {
          sum += s[i] + t[i];
          if (i == len - 1) {
            // There will be ranks left without grids assigned.
            delete [] ret;
            return nullptr;
          } else {
            sum1 = sum + s[i + 1] + t[i + 1];
          }

          if (sum <= avg && avg <= sum1) {
            // Check which is closer to avg
            if (avg - sum <= sum1 - avg) {
              iterator = i;
            } else {
              iterator = i + 1;
              sum = sum1;
            }
            break;
          }
        }

        ret[rank++] = iterator;
        total -= sum;
        iterator++;
      }
    }

    return ret;
  }

  /// Count and sort the local points
  ///
  /// This will assign the local points to the uniform grid, and it will also
  /// rearrange them according to which node of the uniform grid. This will
  /// ultimately return the local counts per uniform grid node which will later
  /// be combined into a global count. The returned counts are allocated in
  /// this routine; the caller assumes ownership of the returned array.
  ///
  /// \param tree - the dual tree
  /// \param p_s - the source data
  /// \param n_sources - the number of sources
  /// \param p_t - the target data
  /// \param n_targets - the number of targets
  /// \param local_offset_s - array that will be allocated and filled with the
  ///                         local offsets for the sources
  /// \param local_offset_t - array that will be allocated and filled with the
  ///                         local offsets for the targets
  ///
  /// \returns - the local counts per uniform grid node
  static int *sort_local_points(DualTree *tree, source_t *p_s,
                                int n_sources, target_t *p_t, int n_targets,
                                int **local_offset_s, int **local_offset_t) {
    int *local_count = new int[tree->dim3_ * 2]();
    int *local_scount = local_count;
    int *local_tcount = &local_count[tree->dim3_];

    int *gid_of_sources = new int[n_sources]();
    int *gid_of_targets = new int[n_targets]();

    // Assign points to the grid
    hpx_addr_t assign_done = hpx_lco_and_new(2);
    assert(assign_done != HPX_NULL);
    DomainGeometry *domarg = &tree->domain_;
    hpx_call(HPX_HERE, sourcetree_t::assign_points_, assign_done,
             &p_s, &n_sources, &domarg, &tree->unif_level_,
             &gid_of_sources, &local_scount);
    hpx_call(HPX_HERE, targettree_t::assign_points_, assign_done,
             &p_t, &n_targets, &domarg, &tree->unif_level_,
             &gid_of_targets, &local_tcount);
    hpx_lco_wait(assign_done);
    hpx_lco_delete_sync(assign_done);

    // Exchange counts
    hpx_lco_set(tree->unif_count_, sizeof(int) * tree->dim3_ * 2, local_count,
                HPX_NULL, HPX_NULL);

    // Group points on the same grid
    hpx_addr_t group_done = hpx_lco_and_new(2);
    assert(group_done != HPX_NULL);
    hpx_call(HPX_HERE, sourcetree_t::group_points_, group_done,
             &p_s, &n_sources, &tree->dim3_, &gid_of_sources, &local_scount,
             &local_offset_s);
    // Only reorder targets if the sources and targets are different
    if (!tree->same_sandt_) {
      hpx_call(HPX_HERE, targettree_t::group_points_, group_done,
               &p_t, &n_targets, &tree->dim3_, &gid_of_targets, &local_tcount,
               &local_offset_t);
    } else {
      hpx_lco_and_set(group_done, HPX_NULL);
    }
    hpx_lco_wait(group_done);
    hpx_lco_delete_sync(group_done);

    delete [] gid_of_sources;
    delete [] gid_of_targets;

    return local_count;
  }

  /// Merge incoming points into the sorted list and spawns adapative paritition
  ///
  /// This action is the 'far side' of the send points message. This will
  /// read the incoming points and merge them with the local data. If the
  /// uniform grid node to which they belong has received all of the points it
  /// is waiting for, this routine will then spawn the adaptive partitioning
  /// of that branch.
  ///
  /// This is a marshalled action, and so the message data is rather opaque.
  ///
  /// \param args - a buffer containing the incoming message.
  /// \parma UNUSED - the size of the message.
  static int recv_points_handler(void *args, size_t UNUSED) {
    hpx_addr_t *rwarg = static_cast<hpx_addr_t *>(args);
    RankWise<dualtree_t> global_tree{*rwarg};
    auto local_tree = global_tree.here();

    // Wait until the buffer is allocated before merging incoming messages
    // We could do this as a call when, but then we need to be aware of the
    // addresses for every rank's LCO. For now, we do this, as it is simpler.
    sourcetree_t *stree = local_tree->source_tree_;
    targettree_t *ttree = local_tree->target_tree_;
    hpx_lco_wait(stree->unif_done_);
    hpx_lco_wait(ttree->unif_done_);

    // TODO: This bit where the message is interpreted might be made easier with
    // ReadBuffer
    int *meta = reinterpret_cast<int *>(static_cast<char *>(args)
                                        + sizeof(hpx_addr_t));
    int first = meta[0];
    int last = meta[1];
    int range = last - first + 1;
    int recv_ns = meta[2];
    int recv_nt = meta[3];
    int *count_s = &meta[4]; // Used only if recv_ns > 0
    int *count_t = count_s + range * (recv_ns > 0); // Used only if recv_nt > 0
    source_t *recv_s =
      reinterpret_cast<source_t *>(static_cast<char *>(args) + sizeof(int) * 4 +
                                sizeof(hpx_addr_t) +
                                sizeof(int) * range * (recv_ns > 0) +
                                sizeof(int) * range * (recv_nt > 0));
    target_t *recv_t = reinterpret_cast<target_t *>(recv_s + recv_ns);

    hpx_addr_t done = hpx_lco_and_new(range * 2);

    if (recv_ns) {
      for (int i = first; i <= last; ++i) {
        int incoming_ns = count_s[i - first];
        if (incoming_ns) {
          sourcenode_t *ns = &stree->unif_grid_[i];
          hpx_call(HPX_HERE, sourcetree_t::merge_points_, done,
                   &recv_s, &ns, &incoming_ns, rwarg);
          recv_s += incoming_ns;
        } else {
          hpx_lco_and_set(done, HPX_NULL);
        }
      }
    } else {
      hpx_lco_and_set_num(done, range, HPX_NULL);
    }

    if (recv_nt) {
      for (int i = first; i <= last; ++i) {
        int incoming_nt = count_t[i - first];
        if (incoming_nt) {
          targetnode_t *nt = &ttree->unif_grid_[i];
          hpx_call(HPX_HERE, targettree_t::merge_points_, done,
                   &recv_t, &nt, &incoming_nt, rwarg);
          recv_t += incoming_nt;
        } else {
          hpx_lco_and_set(done, HPX_NULL);
        }
      }
    } else {
      if (local_tree->same_sandt_) {
        // S == T means do a special version of merge.
        for (int i = first; i <= last; ++i) {
          int incoming_nt = count_s[i - first];
          if (incoming_nt) {
            targetnode_t *nt = &ttree->unif_grid_[i];
            sourcenode_t *ns = &stree->unif_grid_[i];
            hpx_call(HPX_HERE, targettree_t::merge_points_same_s_and_t_,
                     done, &nt, &incoming_nt, &ns, rwarg);
          } else {
            hpx_lco_and_set(done, HPX_NULL);
          }
        }
      } else {
        hpx_lco_and_set_num(done, range, HPX_NULL);
      }
    }

    // Wait until the data has been merged before releasing the parcel
    // NOTE: This 'done' will trigger once points are merged. The action that
    // triggers this will spawn more work, but not in a synchronous way.
    hpx_lco_wait(done);
    hpx_lco_delete_sync(done);
    return HPX_SUCCESS;
  }

  /// Send the points to the remote that will assume ownership
  ///
  /// This action sends points to remote localities that have been assigned
  /// the points during distribution.
  ///
  /// \param rank - the rank to which we are sending
  /// \param count_s - the source counts
  /// \param count_t - that target counts
  /// \param offset_s - the source offsets
  /// \param offset_t - the target offsets
  /// \param sources - the source data
  /// \param targets - the target data
  /// \param rwaddr - the global address of the dual tree
  static int send_points_handler(int rank, int *count_s, int *count_t,
                                 int *offset_s, int *offset_t,
                                 source_t *sources, target_t *targets,
                                 hpx_addr_t rwaddr) {
    RankWise<dualtree_t> global_tree{rwaddr};
    auto local_tree = global_tree.here();

    // Note: all the pointers are local to the calling rank.
    int first = local_tree->first(rank);
    int last = local_tree->last(rank);
    int range = last - first + 1;
    int send_ns = 0;
    int send_nt = 0;
    for (int i = first; i <= last; ++i) {
      send_ns += count_s[i];
      send_nt += count_t[i];
    }

    // Clear out the target sends if S == T
    if (local_tree->same_sandt_) {
      send_nt = 0;
    }

    // Parcel message length
    size_t bytes = sizeof(hpx_addr_t) + sizeof(int) * 4;
    if (send_ns) {
      bytes += sizeof(int) * range + sizeof(source_t) * send_ns;
    }
    if (send_nt) {
      bytes += sizeof(int) * range + sizeof(target_t) * send_nt;
    }

    // Acquire parcel
    hpx_parcel_t *p = hpx_parcel_acquire(nullptr, bytes);
    void *data = hpx_parcel_get_data(p);
    hpx_addr_t *rwarg = static_cast<hpx_addr_t *>(data);
    *rwarg = rwaddr;
    int *meta = reinterpret_cast<int *>(
                              static_cast<char *>(data) + sizeof(hpx_addr_t));
    meta[0] = first;
    meta[1] = last;
    meta[2] = send_ns;
    meta[3] = send_nt;

    int *count = &meta[4];
    if (send_ns) {
      memcpy(count, &count_s[first], sizeof(int) * range);
      count += range;
    }

    if (send_nt) {
      memcpy(count, &count_t[first], sizeof(int) * range);
    }

    char *meta_s = static_cast<char *>(data) + sizeof(hpx_addr_t) +
      sizeof(int) * (4 + range * (send_ns > 0) + range * (send_nt > 0));
    if (send_ns) {
      memcpy(meta_s, &sources[offset_s[first]], sizeof(source_t) * send_ns);
    }

    char *meta_t = meta_s + send_ns * sizeof(source_t);
    if (send_nt) {
      memcpy(meta_t, &targets[offset_t[first]], sizeof(target_t) * send_nt);
    }

    hpx_parcel_set_target(p, HPX_THERE(rank));
    hpx_parcel_set_action(p, recv_points_);
    hpx_parcel_send(p, HPX_NULL);

    return HPX_SUCCESS;
  }

  /// The action responsible for dual tree partitioning
  ///
  /// This action is the target of a broadcast and manages all the work
  /// required to build the distributed trees.
  ///
  /// \param rwtree - the global address of the dual tree
  /// \param sources_gas - the source data
  /// \param targets_gas - the target data
  static int create_dual_tree_handler(hpx_addr_t rwtree,
                                      hpx_addr_t sources_gas,
                                      hpx_addr_t targets_gas) {
    int rank = hpx_get_my_rank();
    int num_ranks = hpx_get_num_ranks();

    RankWise<dualtree_t> global_tree{rwtree};
    auto tree = global_tree.here();

    Array<source_t> sources{sources_gas};
    sourceref_t src_ref = sources.ref(rank);
    sourcearraydata_t src_data = src_ref.pin();
    source_t *p_s = src_data.value();
    int n_sources = src_ref.n();

    Array<target_t> targets{targets_gas};
    targetref_t trg_ref = targets.ref(rank);
    targetarraydata_t trg_data = trg_ref.pin();
    target_t *p_t = trg_data.value();
    int n_targets = trg_ref.n();

    // Assign points to uniform grid
    int *local_offset_s{nullptr};
    int *local_offset_t{nullptr};
    int *local_count = sort_local_points(&*tree, p_s, n_sources,
                                         p_t, n_targets, &local_offset_s,
                                         &local_offset_t);
    int *local_scount = local_count;
    int *local_tcount = &local_count[tree->dim3_];

    // Compute point distribution
    hpx_lco_get(tree->unif_count_, sizeof(int) * (tree->dim3_ * 2),
                tree->unif_count_src());
    tree->distribute_ = distribute_points(num_ranks,
                                          tree->unif_count_src(),
                                          tree->dim3_);

    // Exchange points
    sourcenode_t *ns = tree->source_tree_->unif_grid_;
    targetnode_t *nt = tree->target_tree_->unif_grid_;
    int firstarg = tree->first(rank);
    int lastarg = tree->last(rank);
    sourcetree_t::init_point_exchange(tree->source_tree_, firstarg, lastarg,
        tree->unif_count_src(), local_scount, local_offset_s, &tree->domain_,
        tree->refinement_limit_, p_s, ns);
    if (!tree->same_sandt_) {
      targettree_t::init_point_exchange(tree->target_tree_, firstarg, lastarg,
          tree->unif_count_tar(), local_tcount, local_offset_t, &tree->domain_,
          tree->refinement_limit_, p_t, nt);
    } else {
      targettree_t::init_point_exchange_same_s_and_t(tree->target_tree_,
          firstarg, lastarg, tree->unif_count_tar(), local_tcount,
          local_offset_t, &tree->domain_, tree->refinement_limit_, p_t, nt,
          ns, tree->source_tree_);
    }

    // So this one is pretty simple. It sends those points from this rank
    // going to the other rank in a parcel.
    for (int r = 0; r < num_ranks; ++r) {
      if (r != rank) {
        hpx_call(HPX_HERE, send_points_, HPX_NULL, &r, &local_scount,
                 &local_tcount, &local_offset_s, &local_offset_t,
                 &p_s, &p_t, &rwtree);
      }
    }

    hpx_addr_t dual_tree_complete = hpx_lco_and_new(2 * tree->dim3_);
    assert(dual_tree_complete != HPX_NULL);

    for (int r = 0; r < num_ranks; ++r) {
      int first = tree->first(r);
      int last = tree->last(r);

      if (r == rank) {
        for (int i = first; i <= last; ++i) {
          if (*(tree->unif_count_src(i)) == 0) {
            hpx_lco_and_set(dual_tree_complete, HPX_NULL);
          } else {
            source_t *arg = tree->sorted_src().value();
            int typearg = 0;
            assert(ns[i].complete() != HPX_NULL);
            hpx_call_when_with_continuation(ns[i].complete(),
                HPX_HERE, sourcetree_t::send_node_,
                dual_tree_complete, hpx_lco_set_action,
                &ns, &arg, &i, &rwtree, &typearg);
          }

          if (*(tree->unif_count_tar(i)) == 0) {
            hpx_lco_and_set(dual_tree_complete, HPX_NULL);
          } else {
            target_t *arg = tree->sorted_tar().value();
            int typearg = 1;
            assert(nt[i].complete() != HPX_NULL);
            hpx_call_when_with_continuation(nt[i].complete(),
                HPX_HERE, targettree_t::send_node_,
                dual_tree_complete, hpx_lco_set_action,
                &nt, &arg, &i, &rwtree, &typearg);
          }
        }
      } else {
        for (int i = first; i <= last; ++i) {
          if (*(tree->unif_count_src(i)) == 0) {
            hpx_lco_and_set(dual_tree_complete, HPX_NULL);
          } else {
            assert(ns[i].complete() != HPX_NULL);
            hpx_call_when_with_continuation(ns[i].complete(),
                dual_tree_complete, hpx_lco_set_action,
                ns[i].complete(), hpx_lco_delete_action,
                nullptr, 0);
          }

          if (*(tree->unif_count_tar(i)) == 0) {
            hpx_lco_and_set(dual_tree_complete, HPX_NULL);
          } else {
            assert(nt[i].complete() != HPX_NULL);
            hpx_call_when_with_continuation(nt[i].complete(),
                dual_tree_complete, hpx_lco_set_action,
                nt[i].complete(), hpx_lco_delete_action,
                nullptr, 0);
          }
        }
      }
    }

    hpx_lco_wait(dual_tree_complete);
    hpx_lco_delete_sync(dual_tree_complete);

    // Replace segment in the array
    hpx_addr_t old_src_data = sources.replace(tree->source_tree_->sorted_);
    hpx_gas_free_sync(old_src_data);
    if (!tree->same_sandt_) {
      hpx_addr_t old_tar_data = targets.replace(tree->target_tree_->sorted_);
      hpx_gas_free_sync(old_tar_data);
    }

    delete [] local_count;
    delete [] local_offset_s;
    delete [] local_offset_t;

    return HPX_SUCCESS;
  }

  /// Action to destroy the tree
  ///
  /// This action is the target of a broadcast that is used to destroy the
  /// tree.
  ///
  /// \param rwtree - global address of the dual tree
  static int finalize_partition_handler(hpx_addr_t rwtree) {
    RankWise<dualtree_t> global_tree{rwtree};
    auto tree = global_tree.here();
    tree->clear_data();
    return HPX_SUCCESS;
  }

  //
  // Now for the data members
  //

  DomainGeometry domain_;     /// domain size
  int refinement_limit_;      /// refinement threshold
  int unif_level_;            /// level of uniform partition
  int dim3_;                  /// number of uniform nodes
  hpx_addr_t unif_count_;     /// LCO reducing the uniform counts
  int *unif_count_value_;     /// local data storing the uniform counts
  int *distribute_;           /// the computed distribution of the nodes
  method_t method_;           /// method used during DAG discovery

  sourcetree_t *source_tree_; /// The source tree
  targettree_t *target_tree_; /// The target tree

  int same_sandt_;            /// Made from the same sources and targets


  static hpx_action_t domain_geometry_init_;
  static hpx_action_t domain_geometry_op_;
  static hpx_action_t set_domain_geometry_;
  static hpx_action_t init_partition_;
  static hpx_action_t recv_points_;
  static hpx_action_t send_points_;
  static hpx_action_t create_dual_tree_;
  static hpx_action_t finalize_partition_;
};

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t DualTree<S, T, E, M, D>::domain_geometry_init_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t DualTree<S, T, E, M, D>::domain_geometry_op_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t DualTree<S, T, E, M, D>::set_domain_geometry_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t DualTree<S, T, E, M, D>::init_partition_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t DualTree<S, T, E, M, D>::recv_points_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t DualTree<S, T, E, M, D>::send_points_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t DualTree<S, T, E, M, D>::create_dual_tree_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t DualTree<S, T, E, M, D>::finalize_partition_ = HPX_ACTION_NULL;

#endif
