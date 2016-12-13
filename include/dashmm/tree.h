// =============================================================================
//  This file is part of:
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  DASHMM is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  DASHMM is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with DASHMM. If not, see <http://www.gnu.org/licenses/>.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================


#ifndef __DASHMM_TREE_H__
#define __DASHMM_TREE_H__


/// \file
/// \brief DualTree related types


// C library
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cassert>

// C++ library
#include <algorithm>
#include <functional>
#include <vector>

// HPX-5
#include "hpx/hpx.h"

// DASHMM
#include "dashmm/array.h"
#include "dashmm/dag.h"
#include "dashmm/expansionlco.h"
#include "dashmm/index.h"
#include "dashmm/point.h"
#include "dashmm/domaingeometry.h"
#include "dashmm/rankwise.h"
#include "dashmm/reductionops.h"


namespace dashmm {

// Forward declare Registrars for the objects in this file
template <typename Record>
class NodeRegistrar;

template <typename Source, typename Target, typename Record,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class> class Method>
class TreeRegistrar;

template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class> class Method>
class DualTreeRegistrar;


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
  ///
  /// \param index - the index of the node
  Node(Index index)
      : idx{index}, parts{}, parent{nullptr}, dag{index}, first_{0} {
    for (int i = 0; i < 8; ++i) {
      child[i] = nullptr;
    }
    sema_ = HPX_NULL;
    complete_ = hpx_lco_and_new(8);
  }

  /// Construct with an index, a particle segment and a parent
  ///
  /// This will also create the completion detection LCO.
  ///
  /// \param index - the node index
  /// \param parts - the particles inside the volume represented by this node
  /// \param parent - the parent of this node
  Node(Index index, arrayref_t parts, node_t *parent)
      : idx{index}, parts{parts}, parent{parent}, dag{index}, first_{0} {
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
  /// @p same_sandt will be nonzero only for target nodes, and only sometimes.
  /// In this case, the following does not sort the records, but merely finds
  /// the split point, which will have been established already by the source
  /// tree partioning.
  ///
  /// \param threshold - the partitioning threshold
  /// \param geo - the domain geometry
  /// \param same_sandt - is this a run where the sources and targets are
  ///                     identical.
  void partition(int threshold, DomainGeometry *geo, int same_sandt) {
    size_t num_points = num_parts();
    assert(num_points >= 1);
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
      record_t *p = parts.data();

      // Sort the particles among the children
      record_t *splits[9]{};
      splits[0] = p;
      splits[8] = &p[num_points];

      auto z_comp = [&center_z](record_t &a) {
        return a.position.z() < center_z;
      };
      auto y_comp = [&center_y](record_t &a) {
        return a.position.y() < center_y;
      };
      auto x_comp = [&center_x](record_t &a) {
        return a.position.x() < center_x;
      };

      if (same_sandt) {
        splits[4] = std::partition_point(splits[0], splits[8], z_comp);

        splits[2] = std::partition_point(splits[0], splits[4], y_comp);
        splits[6] = std::partition_point(splits[4], splits[8], y_comp);

        splits[1] = std::partition_point(splits[0], splits[2], x_comp);
        splits[3] = std::partition_point(splits[2], splits[4], x_comp);
        splits[5] = std::partition_point(splits[4], splits[6], x_comp);
        splits[7] = std::partition_point(splits[6], splits[8], x_comp);
      } else {
        splits[4] = std::partition(splits[0], splits[8], z_comp);

        splits[2] = std::partition(splits[0], splits[4], y_comp);
        splits[6] = std::partition(splits[4], splits[8], y_comp);

        splits[1] = std::partition(splits[0], splits[2], x_comp);
        splits[3] = std::partition(splits[2], splits[4], x_comp);
        splits[5] = std::partition(splits[4], splits[6], x_comp);
        splits[7] = std::partition(splits[6], splits[8], x_comp);
      }

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
      curr->dag.set_index(parent->idx.child(which));
      parent->child[which] = curr;
    }
  }

  /// Remove a given child from this node
  void remove_child(const node_t *target) {
    bool found{false};
    for (int i = 0; i < 8; ++i) {
      if (child[i] == target) {
        child[i] = nullptr;
        found = true;
        break;
      }
    }
    assert(found);
  }

  /// This will recurse down to the uniform level and remove pointless nodes
  ///
  /// This returns true if this node has no children after the work is
  /// completed for this node.
  ///
  /// NOTE: This should be a safe recursion because the depth of the top
  /// part of the tree grows very slowly.
  bool remove_downward_links(int limit, int level) {
    if (level < limit) {
      for (int i = 0; i < 8; ++i) {
        if (child[i] && child[i]->remove_downward_links(limit, level + 1)) {
          child[i]->parent = nullptr; // Remove upward link
          child[i] = nullptr;         // Remove downward link
        }
      }
    }
    return (is_leaf() && parts.n() == 0);
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
  DAGInfo dag;                    /// The DAG info for this node

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
                    template <typename, typename> class> class Method>
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
                    template <typename, typename> class> class Method>
class Tree {
 public:
  using record_t = Record;
  using node_t = Node<Record>;
  using sourcenode_t = Node<Source>;
  using targetnode_t = Node<Target>;
  using arrayref_t = ArrayRef<Record>;
  using tree_t = Tree<Source, Target, Record, Expansion, Method>;
  using sourcetree_t = Tree<Source, Target, Source, Expansion, Method>;
  using dualtree_t = DualTree<Source, Target, Expansion, Method>;

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
  ///
  /// \param tree - the address of the tree on which to act
  /// \param unif_level - the uniform partitioning level
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
            scnode->dag.set_index(snode->idx.child(i));
            scnode->parent = snode;
          }
        } else {
          // The final level does need the ordering
          for (int i = 0; i < 8; ++i) {
            Index cindex = snode->idx.child(i);
            uint64_t morton = morton_key(cindex.x(), cindex.y(), cindex.z());
            snode->child[i] = &tree->unif_grid_[morton];
            tree->unif_grid_[morton].idx = cindex;
            tree->unif_grid_[morton].dag.set_index(cindex);
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
  ///
  /// \param tree - the tree on which to act
  /// \param ndim - the size of the uniform grid
  /// \param first - the first owned node of the uniform grid
  /// \param last - the last (inclusive) owned node of the uniform grid
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
  ///
  /// \param P - the records
  /// \param npts - the number of records
  /// \param geo - the overall domain geometry
  /// \param unif_level - the uniform partitioning level
  /// \param gid [out] - the morton key for each record
  /// \param count [out] - the number of points per uniform grid node
  ///
  /// \return HPX_SUCCESS
  static int assign_points_to_unif_grid(const record_t *P, int npts,
                                         const DomainGeometry *geo,
                                         int unif_level, int *gid,
                                         int *count) {
    Point corner = geo->low();
    double scale = 1.0 / geo->size();

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
  ///
  /// \param p_in - the input records; will be sorted
  /// \param npts - the number of records
  /// \param dim3 - the size of the uniform grid
  /// \param gid_of_points - the morton key for the records
  /// \param count - the number of records per bin
  /// \param retval [out] - offsets into the record list for each
  ///
  /// \returns - HPX_SUCCESS
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
  ///
  /// \param message_buffer - the message
  /// \param UNUSED - the size of the buffer
  ///
  /// \returns - HPX_SUCCESS
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
  ///
  /// \param n - the uniform grid nodes
  /// \param sorted - the sorted records
  /// \param id - the uniform grid node in question
  /// \param rwaddr - global address of the DualTree
  /// \param type - indicates source or target tree
  ///
  /// \returns - HPX_SUCCESS
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
    hpx_lco_delete_sync(done);

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

    if (range > 0) {
      // Compute global_offset
      int *global_offset = new int[range]();
      size_t num_points = global_count[first];
      for (int i = first + 1; i <= last; ++i) {
        num_points += global_count[i];
        global_offset[i - first] = global_offset[i - first - 1] +
                                   global_count[i - 1];
      }

      if (num_points > 0) {
        record_t *sorted_data = reinterpret_cast<record_t *>(
                                    new char[num_points * sizeof(record_t)]);
        sorted_ref = arrayref_t{sorted_data, num_points};

        for (int i = first; i <= last; ++i) {
          node_t *curr = &n[i];
          curr->parts = sorted_ref.slice(global_offset[i - first],
                                         global_count[i]);

          if (local_count[i]) {
            // Copy local points before merging remote points
            curr->lock();
            auto sorted = curr->parts.data();
            memcpy(sorted + curr->first(),
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
      delete [] global_offset;
    }

    tree->sorted_ = sorted_ref;
    hpx_lco_and_set(tree->unif_done_, HPX_NULL);
  }

  /// Initialize point exchange between ranks for S == T case
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
  /// \param snodes - the source nodes in the uniform grid
  /// \param source_tree - the source tree
  static void init_point_exchange_same_s_and_t(
      tree_t *tree, int first, int last, const int *global_count,
      const int *local_count, const int *local_offset,
      const DomainGeometry *geo, int threshold, const record_t *temp,
      node_t *n, sourcenode_t *snodes, sourcetree_t *source_tree) {
    int range = last - first + 1;

    if (range > 0) {
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
          curr->parts = arrayref_t{(record_t *)sparts.data(), sparts.n()};

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
      delete [] global_offset;
    }

    // Again, the target tree reuses the data from the source tree
    ArrayRef<Source> ssort = source_tree->sorted();
    tree->sorted_ = arrayref_t{(record_t *)ssort.data(), ssort.n()};
    hpx_lco_and_set(tree->unif_done_, HPX_NULL);
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
  ///
  /// \returns - HPX_SUCCESS
  static int merge_points_handler(record_t *temp, node_t *n, int n_arrived,
                                  hpx_addr_t rwgas) {
    RankWise<dualtree_t> global_tree{rwgas};
    auto local_tree = global_tree.here();

    // Note: all the pointers are local to the calling rank.
    n->lock();
    size_t first = n->first();
    record_t *p = n->parts.data();
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
  ///
  /// \returns - HPX_SUCCESS
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

  /// Compute the uniform grid index for a given Index
  ///
  /// \param idx - the index of the node
  /// \param uniflevel - the uniform level of the tree
  ///
  /// \returns - the index in the uniform grid that gives the ancestor of
  ///            the given Index
  static int get_unif_grid_index(const Index &idx, int uniflevel) {
    int delta = idx.level() - uniflevel;

    if (delta < 0) {
      return -1;
    }

    if (delta > 0) {
      Index tester = idx.parent(delta);
      return morton_key(tester.x(), tester.y(), tester.z());
    } else {
      return morton_key(idx.x(), idx.y(), idx.z());
    }
  }

  /// Find the LCO address for a given index and a given operation
  ///
  /// \param idx - the Index of the node in question
  /// \param op - the edge type connecting to the index in question
  ///
  /// \returns - global address of the LCO serving as target of the edge
  hpx_addr_t lookup_lco_addx(Index idx, Operation op) {
    // This should walk to the node containing the LCO we care about
    node_t *curr = root_;
    bool not_found = true;
    while (not_found) {
      int dlevel = idx.level() - curr->idx.level();
      assert(dlevel >= 0);
      if (dlevel == 0) {
        break;
      }
      assert(curr->idx == idx.parent(dlevel));
      int which = idx.parent(dlevel - 1).which_child();
      curr = curr->child[which];
    }
    assert(curr->idx == idx);

    hpx_addr_t retval{HPX_NULL};
    switch (op) {
      case Operation::Nop:
        assert(0 && "problem in lookup");
        break;
      case Operation::StoM:  // NOTE: fallthrough here
      case Operation::StoL:
      case Operation::MtoM:
      case Operation::MtoL:
      case Operation::LtoL:
        assert(curr->dag.has_normal());
        retval = curr->dag.normal()->global_addx;
        break;
      case Operation::MtoT: // NOTE: fallthrough here
      case Operation::LtoT:
      case Operation::StoT:
        assert(curr->dag.has_parts());
        retval = curr->dag.parts()->global_addx;
        break;
      case Operation::MtoI: // NOTE: fallthrough
      case Operation::ItoI:
        assert(curr->dag.has_interm());
        retval = curr->dag.interm()->global_addx;
        break;
      case Operation::ItoL:
        assert(curr->dag.has_normal());
        retval = curr->dag.normal()->global_addx;
        break;
    }

    return retval;
  }


private:
  friend class DualTree<Source, Target, Expansion, Method>;
  friend class TreeRegistrar<Source, Target, Record, Expansion, Method>;

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
                    template <typename, typename> class> class M>
hpx_action_t Tree<S, T, R, E, M>::setup_basics_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Tree<S, T, R, E, M>::delete_tree_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Tree<S, T, R, E, M>::recv_node_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Tree<S, T, R, E, M>::send_node_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Tree<S, T, R, E, M>::assign_points_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Tree<S, T, R, E, M>::group_points_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Tree<S, T, R, E, M>::merge_points_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Tree<S, T, R, E, M>::merge_points_same_s_and_t_ =
                                                          HPX_ACTION_NULL;


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
                    template <typename, typename> class> class Method>
class DualTree {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, Method>;
  using targetlco_t = TargetLCO<Source, Target, Expansion, Method>;
  using method_t = Method<Source, Target, Expansion>;
  using sourcenode_t = Node<Source>;
  using targetnode_t = Node<Target>;
  using sourceref_t = ArrayRef<Source>;
  using targetref_t = ArrayRef<Target>;
  using sourcetree_t = Tree<Source, Target, Source, Expansion, Method>;
  using targettree_t = Tree<Source, Target, Target, Expansion, Method>;
  using dualtree_t = DualTree<Source, Target, Expansion, Method>;

  /// Construction is always default
  DualTree()
    : domain_{}, refinement_limit_{1}, unif_level_{1}, dim3_{8},
      unif_count_{HPX_NULL}, unif_count_value_{nullptr},
      distribute_{nullptr}, method_{}, source_tree_{nullptr},
      target_tree_{nullptr} { }

  /// We delete the copy constructor and copy assignement operator.
  DualTree(const dualtree_t &other) = delete;
  dualtree_t &operator=(const dualtree_t &other) = delete;

  /// Return the number of source points in the uniform grid
  ///
  /// This routine returns the address of the given information,
  ///
  /// \param i - the uniform grid node in question
  int *unif_count_src(size_t i = 0) const {return &unif_count_value_[i];}

  /// Return the number of target points in the uniform grid
  ///
  /// This routine returns the address of the given information,
  ///
  /// \param i - the uniform grid node in question
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
  size_t sorted_src_count() const {return source_tree_->sorted_.n();}

  /// Return the number of post-sorting targets.
  size_t sorted_tar_count() const {return target_tree_->sorted_.n();}

  /// Return the sorted sources.
  source_t *sorted_src() const {return source_tree_->sorted_.data();}

  /// Return the sorted targets.
  target_t *sorted_tar() const {return target_tree_->sorted_.data();}

  /// Return the method this object will use for DAG operations.
  const method_t &method() const {return method_;}

  /// Set the method this object will use for DAG operations.
  void set_method(const method_t &method) {method_ = method;}

  /// Lookup target LCO address
  ///
  /// \param idx - index of LCO to look up
  /// \param op - operation of the edge leading to the LCO in question
  ///
  /// \returns - global address of LCO
  hpx_addr_t lookup_lco_addx(Index idx, Operation op) {
    bool search_source{true};
    switch(op) {
      case Operation::Nop:
        assert(0 && "Problem in address search");
        break;
      case Operation::StoM:
        break;
      case Operation::StoL:
        search_source = false;
        break;
      case Operation::MtoM:
        break;
      case Operation::MtoL:
        search_source = false;
        break;
      case Operation::LtoL:
        search_source = false;
        break;
      case Operation::MtoT:
        search_source = false;
        break;
      case Operation::LtoT:
        search_source = false;
        break;
      case Operation::StoT:
        search_source = false;
        break;
      case Operation::MtoI:
        break;
      case Operation::ItoI:
        search_source = false;
        break;
      case Operation::ItoL:
        search_source = false;
        break;
    }

    if (search_source) {
      return source_tree_->lookup_lco_addx(idx, op);
    } else {
      return target_tree_->lookup_lco_addx(idx, op);
    }
  }

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
    delete [] rank_map_;
  }

  /// Return the first uniform grid node owned by the given rank.
  int first(int rank) const {return rank == 0 ? 0 : distribute_[rank - 1] + 1;}

  /// Return the last uniform grid node owned by the given rank.
  int last(int rank) const {return distribute_[rank];}

  /// Return the rank owning the given unif grid node
  int rank_of_unif_grid(int idx) const {return rank_map_[idx];}

  /// Create the DAG for this tree using the method specified for this object.
  ///
  /// This will allocate and collect the DAG nodes into the returned object.
  ///
  /// \returns - the resulting DAG.
  DAG *create_DAG() {
    // Do work on the source tree
    hpx_addr_t sdone = hpx_lco_future_new(sizeof(int));
    assert(sdone != HPX_NULL);

    dualtree_t *thetree = this;
    hpx_call(HPX_HERE, source_apply_method_, HPX_NULL,
             &thetree, &source_tree_->root_, &sdone);

    hpx_lco_wait(sdone);
    hpx_lco_delete_sync(sdone);

    // Do work on the target tree
    hpx_addr_t tdone = hpx_lco_and_new(1);
    assert(tdone != HPX_NULL);

    std::vector<sourcenode_t *> *consider = new std::vector<sourcenode_t *>{};
    consider->push_back(source_tree_->root_);
    targetnode_t *trgaddx = target_tree_->root_;
    hpx_call(HPX_HERE, target_apply_method_, HPX_NULL,
             &thetree, &trgaddx, &consider, &same_sandt_, &tdone);

    hpx_lco_wait(tdone);
    hpx_lco_delete_sync(tdone);

    DAG *retval = new DAG{};
    collect_DAG_nodes(retval);
    return retval;
  }

  /// Traverse the tree and collect the DAG nodes
  ///
  /// This will collect the nodes of the DAG into three groups: the source
  /// nodes of the DAG, the target nodes of the DAG, and all other nodes.
  /// The source and target nodes are the input and output nodes respectively.
  /// The remainder are those nodes containing an intermediate computation.
  ///
  /// This is a synchronous operation.
  ///
  /// \param dag - a DAG object to be populated
  void collect_DAG_nodes(DAG *dag) {
    collect_DAG_nodes_from_S_node(source_tree_->root_, dag->source_leaves,
                                  dag->source_nodes);
    collect_DAG_nodes_from_T_node(target_tree_->root_, dag->target_leaves,
                                  dag->target_nodes);
  }

  /// Create the LCOs from the DAG
  ///
  /// This will traverse the source and target tree creating any needed
  /// expansion LCOs. Further, it will create the target LCOs in the target
  /// tree.
  ///
  /// This is a synchronous operation.
  void create_expansions_from_DAG(hpx_addr_t rwtree) {
    hpx_addr_t done = hpx_lco_and_new(2);
    assert(done != HPX_NULL);

    dualtree_t *argthis = this;
    sourcenode_t *srcaddx = source_tree_->root_;
    hpx_call(HPX_HERE, create_S_expansions_from_DAG_, HPX_NULL,
             &done, &argthis, &srcaddx, &rwtree);
    targetnode_t *trgaddx = target_tree_->root_;
    hpx_call(HPX_HERE, create_T_expansions_from_DAG_, HPX_NULL,
             &done, &argthis, &trgaddx, &rwtree);

    hpx_lco_wait(done);
    hpx_lco_delete_sync(done);
  }

  /// Sets up termination detection for a DASHMM evaluation
  ///
  /// This is an asynchronous operation. The returned LCO becomes the
  /// responsibility of the caller.
  ///
  /// \param targets - vector of target nodes in the DAG
  /// \param internals - vector of internal nodes in the DAG
  ///
  /// \returns - LCO that will signal that all targets have completed their
  ///            computation.
  hpx_addr_t setup_termination_detection(DAG *dag) {
    size_t n_targs = dag->target_leaves.size();
    size_t n_tinternal = dag->target_nodes.size();
    size_t n_sinternal = dag->source_nodes.size();

    hpx_addr_t retval = hpx_lco_and_new(n_targs + n_tinternal + n_sinternal);
    assert(retval != HPX_NULL);

    std::vector<DAGNode *> *argaddx = &dag->target_leaves;
    std::vector<DAGNode *> *tinternalsaddx = &dag->target_nodes;
    std::vector<DAGNode *> *sinternalsaddx = &dag->source_nodes;
    hpx_call(HPX_HERE, termination_detection_, HPX_NULL, &retval,
             &argaddx, &n_targs, &tinternalsaddx, &n_tinternal,
             &sinternalsaddx, &n_sinternal);

    return retval;
  }

  /// Sets the edge lists of the expansion LCOs
  ///
  /// This is a separate phase as the DAG is constructed before the LCOs
  /// serving the work for a DAG node are created. Once the expansion LCOs
  /// are created, the addresses can then be shared with any LCO that needs
  /// that information.
  ///
  /// This is an asynchronous operation. The work will have started when this
  /// function returns, but it may not have ended. There is no returned LCO
  /// as the termination detection cannot trigger before this is done.
  ///
  /// \param dag - DAG object
  void setup_edge_lists(DAG *dag) {
    DAGNode **sdata = dag->source_nodes.data();
    size_t n_snodes = dag->source_nodes.size();
    DAGNode **tdata = dag->target_nodes.data();
    size_t n_tnodes = dag->target_nodes.size();
    hpx_call(HPX_HERE, edge_lists_, HPX_NULL,
             &sdata, &n_snodes, &tdata, &n_tnodes);
  }

  /// Initiate the DAG evaluation
  ///
  /// This starts the work of the evalution by starting the S->* work at the
  /// source nodes of the DAG.
  ///
  /// This is an asynchronous operation. The termination detection cannot
  /// possibly trigger before this is completed, so waiting on the termination
  /// of the full evaluation implicitly waits on this operation.
  ///
  /// \param global_tree - the Dual Tree
  void start_DAG_evaluation(RankWise<dualtree_t> &global_tree) {
    hpx_addr_t rwaddr = global_tree.data();
    hpx_call(HPX_HERE, instigate_dag_eval_, HPX_NULL,
             &rwaddr, &source_tree_->root_);
  }

  /// Destroys the LCOs associated with the DAG
  ///
  /// This is a synchronous operation. This destroys not only the expansion
  /// LCOs, but also the target LCOs.
  ///
  /// \param dag - the DAG
  void destroy_DAG_LCOs(DAG &dag) {
    hpx_addr_t done = hpx_lco_and_new(3);
    assert(done != HPX_NULL);

    DAGNode **data = dag.target_leaves.data();
    size_t n_data = dag.target_leaves.size();
    int type = 1;
    hpx_call(HPX_HERE, destroy_DAG_LCOs_, done, &data, &n_data, &type);

    data = dag.target_nodes.data();
    n_data = dag.target_nodes.size();
    type = 0;
    hpx_call(HPX_HERE, destroy_DAG_LCOs_, done, &data, &n_data, &type);

    data = dag.source_nodes.data();
    n_data = dag.source_nodes.size();
    type = 0;
    hpx_call(HPX_HERE, destroy_DAG_LCOs_, done, &data, &n_data, &type);

    hpx_lco_wait(done);
    hpx_lco_delete_sync(done);
  }


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
  static void destroy(RankWise<dualtree_t> &global_tree) {
    hpx_addr_t rwtree = global_tree.data();
    hpx_bcast_rsync(finalize_partition_, &rwtree);

    auto tree = global_tree.here();
    hpx_lco_delete_sync(tree->unif_count_);

    global_tree.destroy();
  }


 private:
  friend class DualTreeRegistrar<Source, Target, Expansion, Method>;

  /// Edge record for DAG instigation
  struct DAGInstigationRecord {
    Operation op;
    hpx_addr_t target;
    size_t n_parts;
    Index idx;
  };

  /// Action to set the domain geometry given the sources and targets
  ///
  /// This action is the target of a broadcast, and computes the domain for
  /// the local sources and targets. The result is then given to a reduction
  /// LCO which reduces each rank's portion.
  ///
  /// \param sources_gas - the global address of the source data
  /// \param targets_gas - the global address of the target data
  /// \param domain_geometry - a reduction LCO to which the local domain is sent
  ///
  /// \returns HPX_SUCCESS
  static int set_domain_geometry_handler(hpx_addr_t sources_gas,
                                         hpx_addr_t targets_gas,
                                         hpx_addr_t domain_geometry,
                                         int same_sandt) {
    Array<source_t> sources{sources_gas};
    sourceref_t src_ref = sources.ref();
    Source *s = src_ref.data();

    double var[6] = {1e50, -1e50, 1e50, -1e50, 1e50, -1e50};

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
      Target *t = trg_ref.data();

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
  ///
  /// \param values - the LCO's data
  /// \param UNUSED - the size of the data
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
  ///
  /// \param lhs - the LCO data
  /// \param rhs - the set data
  /// \param UNUSED - the size of the set data
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
  /// \param same_sandt - is S == T for this tree
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
  /// \param same_sandt - is S == T for this tree
  ///
  /// \returns - HPX_SUCCESS
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
      new Tree<Source, Target, Source, Expansion, Method>{};
    tree->target_tree_ =
      new Tree<Source, Target, Target, Expansion, Method>{};
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
  /// \param same_sandt - is S == T for this tree
  ///
  /// \returns - the Dual Tree
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
  ///
  /// \returns - the distribution
  static int *distribute_points(int num_ranks, const int *global, int len) {
    int *ret = new int[num_ranks]();

    const int *s = global; // Source counts
    const int *t = &global[len]; // Target counts

    int *cumulative = new int[len + 1];
    cumulative[0] = 0;
    for (int i = 1; i <= len; ++i) {
      cumulative[i] = s[i-1] + t[i-1] + cumulative[i - 1];
    }

    int target_share = cumulative[len] / num_ranks;

    for (int split = 1; split < num_ranks; ++split) {
      int my_target = target_share * split;
      auto fcn = [&my_target] (const int &a) -> bool {return a < my_target;};
      int *splitter = std::partition_point(cumulative, &cumulative[len + 1],
                                           fcn);
      int upper_delta = *splitter - my_target;
      assert(upper_delta >= 0);
      int lower_delta = my_target - *(splitter - 1);
      assert(lower_delta >= 0);
      int split_index = upper_delta < lower_delta
                          ? splitter - cumulative
                          : (splitter - cumulative) - 1;
      --split_index;
      assert(split_index >= 0);

      ret[split - 1] = split_index;
    }
    ret[num_ranks - 1] = len - 1;

    delete [] cumulative;

    return ret;
  }

  /// Generate the mapping from uniform grid node to locality
  ///
  /// \param num_ranks - the number of localities
  void generate_rank_map(int num_ranks) {
    rank_map_ = new int [dim3_];
    assert(rank_map_);

    for (int i = 0; i < num_ranks; ++i) {
      int b = first(i);
      int e = last(i);
      for (int j = b; j <= e; ++j) {
        rank_map_[j] = i;
      }
    }
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
      if (range) {
        hpx_lco_and_set_num(done, range, HPX_NULL);
      }
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
      if (local_tree->same_sandt_ && recv_ns) {
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
        if (range) {
          hpx_lco_and_set_num(done, range, HPX_NULL);
        }
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

  /// This will prune links in the top part of the tree that are not needed
  ///
  /// After exchanging counts, there could be some nodes at the uniform level
  /// that will not have any particles under them. In this case, the links to
  /// these nodes from higher levels should be culled to avoid sending Methods
  /// into pointless work, and potentially scheduling pointless work.
  ///
  /// TODO: At the moment, this is done in serial. It might be useful to
  /// get some parallelism here. Probably the easiest is to do this work
  /// as its own task.
  void prune_topnodes() {
    // Go over uniform level and remove link to any no particle nodes
    int *s_counts = unif_count_src();
    int *t_counts = unif_count_tar();
    sourcenode_t *s_nodes = source_tree_->unif_grid_;
    targetnode_t *t_nodes = target_tree_->unif_grid_;
    for (int i = 0; i < dim3_; ++i) {
      if (s_counts[i] == 0) {
        s_nodes[i].parent->remove_child(&s_nodes[i]);
      }
      if (t_counts[i] == 0) {
        t_nodes[i].parent->remove_child(&t_nodes[i]);
      }
    }

    // Now descend from above and remove any empty nodes, stopping at
    // unif_level_ - 1
    source_tree_->root_->remove_downward_links(unif_level_ - 1, 0);
    target_tree_->root_->remove_downward_links(unif_level_ - 1, 0);
  }

  /// The action responsible for dual tree partitioning
  ///
  /// This action is the target of a broadcast and manages all the work
  /// required to build the distributed trees.
  ///
  /// \param rwtree - the global address of the dual tree
  /// \param sources_gas - the source data
  /// \param targets_gas - the target data
  ///
  /// \return HPX_SUCCESS
  static int create_dual_tree_handler(hpx_addr_t rwtree,
                                      hpx_addr_t sources_gas,
                                      hpx_addr_t targets_gas) {
    int rank = hpx_get_my_rank();
    int num_ranks = hpx_get_num_ranks();

    RankWise<dualtree_t> global_tree{rwtree};
    auto tree = global_tree.here();

    Array<source_t> sources{sources_gas};
    Array<target_t> targets{targets_gas};

    {
      sourceref_t src_ref = sources.ref();
      source_t *p_s = src_ref.data();
      int n_sources = src_ref.n();

      targetref_t trg_ref = targets.ref();
      target_t *p_t = trg_ref.data();
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
      tree->generate_rank_map(num_ranks);

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
            tree->unif_count_tar(), local_tcount, local_offset_t,
            &tree->domain_, tree->refinement_limit_, p_t, nt);
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
              source_t *arg = tree->sorted_src();
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
              target_t *arg = tree->sorted_tar();
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

      // This will prune pointless nodes from the top of the tree.
      tree->prune_topnodes();

      delete [] local_count;
      delete [] local_offset_s;
      delete [] local_offset_t;
    }

    // Replace segment in the array
    char *old_src_data = sources.replace(tree->source_tree_->sorted_);
    if (old_src_data != nullptr) {
      delete [] old_src_data;
    }
    if (!tree->same_sandt_) {
      char *old_tar_data = targets.replace(tree->target_tree_->sorted_);
      if (old_tar_data != nullptr) {
        delete [] old_tar_data;
      }
    }

    return HPX_SUCCESS;
  }

  /// Action to destroy the tree
  ///
  /// This action is the target of a broadcast that is used to destroy the
  /// tree.
  ///
  /// \param rwtree - global address of the dual tree
  ///
  /// \returns - HPX_SUCCESS
  static int finalize_partition_handler(hpx_addr_t rwtree) {
    RankWise<dualtree_t> global_tree{rwtree};
    auto tree = global_tree.here();
    tree->clear_data();
    return HPX_SUCCESS;
  }

  /// Collect DAG nodes from source Tree Nodes
  ///
  /// \param root - tree node
  /// \param sources [out] - DAG nodes that are sources
  /// \param internals [out] - other source tree DAG nodes
  void collect_DAG_nodes_from_S_node(sourcenode_t *root,
                                     std::vector<DAGNode *> &sources,
                                     std::vector<DAGNode *> &internals) {
    for (int i = 0; i < 8; ++i) {
      if (root->child[i]) {
        collect_DAG_nodes_from_S_node(root->child[i], sources, internals);
      }
    }
    root->dag.collect_DAG_nodes(sources, internals);
  }

  /// Collect DAG nodes from target Tree Nodes
  ///
  /// \param root - tree node
  /// \param targets [out] - DAG nodes that are targets
  /// \param internals [out] - other target tree DAG nodes
  void collect_DAG_nodes_from_T_node(targetnode_t *root,
                                     std::vector<DAGNode *> &targets,
                                     std::vector<DAGNode *> &internals) {
    for (int i = 0; i < 8; ++i) {
      if (root->child[i]) {
        collect_DAG_nodes_from_T_node(root->child[i], targets, internals);
      }
    }
    root->dag.collect_DAG_nodes(targets, internals);
  }

  /// Action to create LCOs from the DAG
  ///
  /// This will walk through the source tree and create the Expansion LCOs
  /// needed by @p node. This will spawns recursively as HPX-5 actions the
  /// full tree.
  ///
  /// \param done - LCO address for indicating completion of LCO creation
  /// \param tree - the DualTree
  /// \param node - the node under examination
  /// \param rwtree - the global address of the DualTree
  ///
  /// \returns - HPX_SUCCESS
  static int create_S_expansions_from_DAG_handler(hpx_addr_t done,
                                                  dualtree_t *tree,
                                                  sourcenode_t *node,
                                                  hpx_addr_t rwtree) {
    Point n_center = tree->domain_.center_from_index(node->idx);

    int myrank = hpx_get_my_rank();

    // create the normal expansion if needed
    if (node->dag.has_normal() && node->dag.normal()->locality == myrank) {
      std::unique_ptr<expansion_t> input_expand{
        new expansion_t{n_center, expansion_t::compute_scale(node->idx),
                        kSourcePrimary}
      };
      expansionlco_t expand(node->dag.normal()->in_edges.size(),
                            node->dag.normal()->out_edges.size(),
                            tree->domain_, node->idx, std::move(input_expand),
                            rwtree);
      node->dag.set_normal_expansion(expand.lco());
    }

    // If there is to be an intermediate expansion, create that
    if (node->dag.has_interm() && node->dag.interm()->locality == myrank) {
      std::unique_ptr<expansion_t> interm_expand{
        new expansion_t{n_center, expansion_t::compute_scale(node->idx),
                        kSourceIntermediate}
      };
      expansionlco_t intexp_lco(node->dag.interm()->in_edges.size(),
                                node->dag.interm()->out_edges.size(),
                                tree->domain_, node->idx,
                                std::move(interm_expand),
                                rwtree);
      node->dag.set_interm_expansion(intexp_lco.lco());
    }

    // spawn work at children
    int n_children{node->n_children()};

    if (n_children) {
      hpx_addr_t cdone = hpx_lco_and_new(n_children);
      assert(cdone != HPX_NULL);

      for (int i = 0; i < 8; ++i) {
        if (node->child[i] != nullptr) {
          sourcenode_t *srcaddx = node->child[i];
          hpx_call(HPX_HERE, create_S_expansions_from_DAG_, HPX_NULL,
                   &cdone, &tree, &srcaddx, &rwtree);
        }
      }

      // This will set the parent's LCO as well as delete cdone
      hpx_call_when(cdone, cdone, hpx_lco_delete_action,
                    done, nullptr, 0);
    } else {
      hpx_lco_set(done, 0, nullptr, HPX_NULL, HPX_NULL);
    }

    return HPX_SUCCESS;
  }

  /// Action to create LCOs from the DAG
  ///
  /// This will walk through the target tree and create the Expansion LCOs
  /// needed by @p node. This will spawns recursively as HPX-5 actions the
  /// full tree.
  ///
  /// \param done - LCO address for indicating completion of LCO creation
  /// \param tree - the DualTree
  /// \param node - the node under examination
  /// \param rwtree - the global address of the DualTree
  ///
  /// \returns - HPX_SUCCESS
  static int create_T_expansions_from_DAG_handler(hpx_addr_t done,
                                                  dualtree_t *tree,
                                                  targetnode_t *node,
                                                  hpx_addr_t rwtree) {
    Point n_center = tree->domain_.center_from_index(node->idx);

    int myrank = hpx_get_my_rank();

    // create the normal expansion if needed
    if (node->dag.has_normal() && node->dag.normal()->locality == myrank) {
      std::unique_ptr<expansion_t> input_expand{
        new expansion_t{n_center, expansion_t::compute_scale(node->idx),
                        kTargetPrimary}
      };
      expansionlco_t expand(node->dag.normal()->in_edges.size(),
                            node->dag.normal()->out_edges.size(),
                            tree->domain_, node->idx, std::move(input_expand),
                            rwtree);
      node->dag.set_normal_expansion(expand.lco());
    }

    // If there is to be an intermediate expansion, create that
    if (node->dag.has_interm() && node->dag.interm()->locality == myrank) {
      std::unique_ptr<expansion_t> interm_expand{
        new expansion_t{n_center, expansion_t::compute_scale(node->idx),
                        kTargetIntermediate}
      };
      expansionlco_t intexp_lco(node->dag.interm()->in_edges.size(),
                                node->dag.interm()->out_edges.size(),
                                tree->domain_, node->idx,
                                std::move(interm_expand),
                                rwtree);
      node->dag.set_interm_expansion(intexp_lco.lco());
    }

    // NOTE: this spawn through the tree does not end when the tree ends.
    // Instead, we have to check if this node has a parts node in the DAG.
    // If so, this branch is done, and we need not spawn more.

    // Here is where we make the target lco if needed
    if (node->dag.has_parts()) {
      if (node->dag.parts()->locality == myrank) {
        targetlco_t tlco{node->dag.parts()->in_edges.size(), node->parts};
        node->dag.set_targetlco(tlco.lco(), tlco.n());
      }

      hpx_lco_set(done, 0, nullptr, HPX_NULL, HPX_NULL);
    } else {
      hpx_addr_t cdone = hpx_lco_and_new(node->n_children());
      assert(cdone != HPX_NULL);

      for (int i = 0; i < 8; ++i) {
        if (node->child[i] != nullptr) {
          targetnode_t *trgaddx = node->child[i];
          hpx_call(HPX_HERE, create_T_expansions_from_DAG_, HPX_NULL,
                   &cdone, &tree, &trgaddx, &rwtree);
        }
      }

      // This will set the parent's LCO as well as delete cdone
      hpx_call_when(cdone, cdone, hpx_lco_delete_action,
                    done, nullptr, 0);
    }

    return HPX_SUCCESS;
  }

  /// Action to set the edge lists of the LCOs
  ///
  /// \param snodes - source DAG nodes
  /// \param n_snodes - the number of source nodes
  /// \param tnodes - target DAG nodes
  /// \param n_tnodes - the number of target nodes
  ///
  /// \returns - HPX_SUCCESS
  static int edge_lists_handler(DAGNode **snodes, size_t n_snodes,
                                DAGNode **tnodes, size_t n_tnodes) {
    int myrank = hpx_get_my_rank();
    for (size_t i = 0; i < n_snodes; ++i) {
      if (snodes[i]->locality == myrank) {
        expansionlco_t expand{snodes[i]->global_addx};
        expand.set_out_edge_data(snodes[i]->out_edges);
      }
    }
    for (size_t i = 0; i < n_tnodes; ++i) {
      if (tnodes[i]->locality == myrank) {
        expansionlco_t expand{tnodes[i]->global_addx};
        expand.set_out_edge_data(tnodes[i]->out_edges);
      }
    }
    return HPX_SUCCESS;
  }

  /// Action to start DAG evaluation work
  ///
  /// This is a recursove spawn through the source tree to begin the work of
  /// the DAG evaluation. This will cause the source DAG nodes to begin their
  /// out edges.
  ///
  /// \param rwtree - the global address of the DualTree
  /// \param node - the node being worked on for this action
  ///
  /// \returns - HPX_SUCCESS
  static int instigate_dag_eval_handler(hpx_addr_t rwtree, sourcenode_t *node) {
    RankWise<dualtree_t> global_tree{rwtree};
    auto tree = global_tree.here();

    DAGNode *parts{nullptr};

    int n_children{node->n_children()};
    if (n_children > 0) {
      for (int i = 0; i < 8; ++i) {
        if (node->child[i] != nullptr) {
          hpx_call(HPX_HERE, instigate_dag_eval_, HPX_NULL,
                  &rwtree, &node->child[i]);
        }
      }
    } else {
      // At a leaf, we do actual work
      parts = node->dag.parts();
      if (parts != nullptr && parts->locality != hpx_get_my_rank()) {
        parts = nullptr;
      }
    }

    if (parts) {
      sourceref_t sources = node->parts;

      // We first sort the out edges by locality
      std::sort(parts->out_edges.begin(), parts->out_edges.end(),
                DAG::compare_edge_locality);

      // Make scratch space for the sends
      size_t source_size = sizeof(Source) * sources.n();
      size_t header_size = source_size + sizeof(size_t)
          + sizeof(hpx_addr_t);
      size_t total_size = header_size + sizeof(size_t)
          + parts->out_edges.size() * sizeof(DAGInstigationRecord);
      char *scratch = new char [total_size];
      assert(scratch != nullptr);

      // Copy source data
      auto sref = sources.data();
      {
        WriteBuffer headdata(scratch, header_size);
        assert(headdata.write(sources.n()));

        ReadBuffer sourcedata((char *)sref, source_size);
        assert(headdata.write(sourcedata));

        assert(headdata.write(rwtree));
      }

      int my_rank = hpx_get_my_rank();
      auto begin = parts->out_edges.begin();
      auto end = parts->out_edges.end();
      while (begin != end) {
        int curr_rank = begin->target->locality;
        auto curr = begin;
        while (curr != end && curr->target->locality == curr_rank) {
          ++curr;
        }

        //copy in edge data
        char *edgedata = scratch + header_size;
        size_t *edgecount = reinterpret_cast<size_t *>(edgedata);
        *edgecount = curr - begin;

        DAGInstigationRecord *edgerecords
            = reinterpret_cast<DAGInstigationRecord *>(edgedata
                                                        + sizeof(size_t));
        int i = 0;
        for (auto loop = begin; loop != curr; ++loop) {
          edgerecords[i].op = loop->op;
          edgerecords[i].target = loop->target->global_addx;
          edgerecords[i].n_parts = loop->target->n_parts;
          edgerecords[i].idx = loop->target->idx;
          ++i;
        }

        // Send parcel or do the work
        if (curr_rank == my_rank) {
          instigate_dag_eval_work(sources.n(), sref, tree->domain_,
                                  *edgecount, edgerecords);
        } else {
          size_t parcel_size = header_size + sizeof(size_t)
                               + sizeof(DAGInstigationRecord) * (*edgecount);
          hpx_parcel_t *parc = hpx_parcel_acquire(scratch, parcel_size);
          hpx_parcel_set_action(parc, instigate_dag_eval_remote_);
          hpx_parcel_set_target(parc, HPX_THERE(curr_rank));

          hpx_parcel_send_sync(parc);
        }

        begin = curr;
      }

      delete [] scratch;
    }

    return HPX_SUCCESS;
  }

  /// Action on remote side for DAG instigation
  ///
  /// Similar to the out edges for the Expansion LCOs, the edges out of the
  /// source nodes are sorted and the data is sent across the network once to
  /// each locality. This action handles the fan out once the data reaches
  /// the target locality.
  ///
  /// \param message - the message data
  /// \param bytes - the message size
  ///
  /// \returns - HPX_SUCCESS
  static int instigate_dag_eval_remote_handler(char *message, size_t bytes) {
    // unpack message into arguments to the local work function
    auto input = ReadBuffer(message, bytes);
    size_t n_src{};
    input.read(&n_src);

    Source *sources = input.interpret_array<Source>(n_src);

    hpx_addr_t tree_addx{};
    input.read(&tree_addx);
    RankWise<dualtree_t> global_tree{tree_addx};
    auto local_tree = global_tree.here();

    size_t n_edges{};
    input.read(&n_edges);
    DAGInstigationRecord *edges
        = input.interpret_array<DAGInstigationRecord>(n_edges);

    // Detect if the edges have unknown target addresses and lookup the
    // correct address
    for (size_t i = 0; i < n_edges; ++i) {
      if (edges[i].target == HPX_NULL) {
        edges[i].target = local_tree->target_tree_->lookup_lco_addx(
                                edges[i].idx, edges[i].op);
      }
    }

    instigate_dag_eval_work(n_src, sources, local_tree->domain_,
                            n_edges, edges);

    return HPX_SUCCESS;
  }

  /// Perform the actual work of DAG instigation
  ///
  /// \param n_src - the number of sources
  /// \param sources - the source records
  /// \param domain - the domain geometry
  /// \param n_edges - the number of edges to process
  /// \param edge - the edge data
  static void instigate_dag_eval_work(size_t n_src, Source *sources,
                                      DomainGeometry &domain,
                                      size_t n_edges,
                                      DAGInstigationRecord *edge) {
    // loop over edges
    for (size_t i = 0; i < n_edges; ++i) {
      switch (edge[i].op) {
        case Operation::Nop:
          assert(0 && "Trouble handling DAG instigation");
          break;
        case Operation::StoM:
          {
            expansionlco_t expand{edge[i].target};
            Point center = domain.center_from_index(edge[i].idx);
            expand.S_to_M(center, sources, n_src, edge[i].idx);
          }
          break;
        case Operation::StoL:
          {
            expansionlco_t expand{edge[i].target};
            Point center = domain.center_from_index(edge[i].idx);
            expand.S_to_L(center, sources, n_src, edge[i].idx);
          }
          break;
        case Operation::MtoM:   // NOTE: Fall-through
        case Operation::MtoL:   //   |
        case Operation::LtoL:   //   |
        case Operation::MtoT:   //   |
        case Operation::LtoT:   //   v
          assert(0 && "Trouble handling DAG instigation");
          break;
        case Operation::StoT:
          {
            // S_to_T on expansion LCOs do not need any of the
            // expansionlco_t's state, so we create a default object.
            expansionlco_t expand{HPX_NULL};
            targetlco_t targets{edge[i].target, edge[i].n_parts};
            expand.S_to_T(sources, n_src, targets);
          }
          break;
        default:
          assert(0 && "Trouble handling DAG instigation");
          break;
      }
    }
  }

  /// Action to apply Method::aggregate
  ///
  /// This is an action that is called dependent on the children of this
  /// node all completing their Method application work.
  ///
  /// \param tree - the DualTree
  /// \param node - the node of the source tree
  /// \param cdone - the completion detection LCO that is deleted by this action
  /// \param done - the completion LCO that is set by this action
  ///
  /// \returns - HPX_SUCCESS
  static int source_apply_method_child_done_handler(dualtree_t *tree,
                                                    sourcenode_t *node,
                                                    hpx_addr_t cdone,
                                                    hpx_addr_t done) {
    tree->method_.aggregate(node, &tree->domain_);
    int loc{0};
    if (node->idx.level() >= tree->unif_level_) {
      int dag_idx = sourcetree_t::get_unif_grid_index(node->idx,
                                                      tree->unif_level_);
      loc = tree->rank_of_unif_grid(dag_idx);
    }

    int height;
    hpx_lco_get(cdone, sizeof(int), &height);
    hpx_lco_set(done, sizeof(int), &height, HPX_NULL, HPX_NULL);
    hpx_lco_delete_sync(cdone);

    method_t::distropolicy_t::assign_for_source(node->dag, loc, height);
    return HPX_SUCCESS;
  }

  /// Action to apply Method::generate
  ///
  /// This is a recursive spawn through the source tree.
  ///
  /// \param tree - the DualTree
  /// \param node - the node of the source tree
  /// \param done - LCO to signal when application of the method is complete
  ///
  /// \returns - HPX_SUCCESS
  static int source_apply_method_handler(dualtree_t *tree, sourcenode_t *node,
                                         hpx_addr_t done) {
    int n_children = node->n_children();
    if (n_children == 0) {
      tree->method_.generate(node, &tree->domain_);

      int dag_idx = sourcetree_t::get_unif_grid_index(node->idx,
                                                      tree->unif_level_);
      assert(dag_idx >= 0);
      int dag_rank = tree->rank_of_unif_grid(dag_idx);
      node->dag.set_parts_locality(dag_rank);
      node->dag.set_normal_locality(dag_rank);

      method_t::distropolicy_t::assign_for_source(node->dag, dag_rank, 0);

      int height = 0;
      hpx_lco_set(done, sizeof(int), &height, HPX_NULL, HPX_NULL);

      return HPX_SUCCESS;
    }

    hpx_addr_t cdone = hpx_lco_reduce_new(n_children, sizeof(int),
                                          int_max_ident_op, int_max_op);
    assert(cdone != HPX_NULL);

    for (int i = 0; i < 8; ++i) {
      if (node->child[i] == nullptr) continue;
      hpx_call(HPX_HERE, source_apply_method_, HPX_NULL,
               &tree, &node->child[i], &cdone);
    }

    // Once the children are done, call aggregate here, continuing a set to
    // done once that has happened.
    assert(cdone != HPX_NULL);
    hpx_call_when(cdone, HPX_HERE, source_apply_method_child_done_, HPX_NULL,
                  &tree, &node, &cdone, &done);

    return HPX_SUCCESS;
  }

  /// Action to apply Method::inherit and Method::process
  ///
  /// This is a parallel spawn through the tree to apply the method to the
  /// given tree. This will generate the DAG.
  ///
  /// \param tree - the DualTree
  /// \param node - the target node in question
  /// \param consider - the list of source nodes to consider in process
  /// \param same_sandt - is S == T for this evaluation
  /// \param done - LCO to signal with completion
  ///
  /// \returns - HPX_SUCCESS
  static int target_apply_method_handler(dualtree_t *tree, targetnode_t *node,
                                         std::vector<sourcenode_t *> *consider,
                                         int same_sandt, hpx_addr_t done) {
    bool refine = false;
    if (node->idx.level() < tree->unif_level_) {
      refine = true;
    } else if (node->n_children()) {
      refine = tree->method_.refine_test((bool)same_sandt, node, *consider);
    }

    tree->method_.inherit(node, &tree->domain_, !refine);
    tree->method_.process(node, *consider, !refine, &tree->domain_);

    if (!refine) {
      // If we are not refining, then we are at a target leaf. This means
      // we should set the particle and normal DAG nodes to have this locality.

      int dag_idx = targettree_t::get_unif_grid_index(node->idx,
                                                      tree->unif_level_);
      int dag_rank = tree->rank_of_unif_grid(dag_idx);
      node->dag.set_parts_locality(dag_rank);
      node->dag.set_normal_locality(dag_rank);

      hpx_lco_set_lsync(done, 0, nullptr, HPX_NULL);
    } else {
      int n_children = node->n_children();
      hpx_addr_t cdone = hpx_lco_and_new(n_children);
      assert(cdone);

      for (int i = 0; i < 8; ++i) {
        if (node->child[i] == nullptr) continue;
        std::vector<sourcenode_t *> *ccons =
            new std::vector<sourcenode_t *>{};
        *ccons = *consider;
        hpx_call(HPX_HERE, target_apply_method_, HPX_NULL,
                 &tree, &node->child[i], &ccons, &same_sandt, &cdone);
      }

      int loc{0};
      if (node->idx.level() >= tree->unif_level_) {
        int dag_idx = targettree_t::get_unif_grid_index(node->idx,
                                                        tree->unif_level_);
        loc = tree->rank_of_unif_grid(dag_idx);
      }
      method_t::distropolicy_t::assign_for_target(node->dag, loc);

      assert(cdone != HPX_NULL);
      hpx_call_when(cdone, cdone, hpx_lco_delete_action,
                    done, nullptr, 0);
    }

    delete consider;

    return HPX_SUCCESS;
  }

  /// Action to set up termination detection for the DAG evaluation
  ///
  /// \param done - LCO for termination detection
  /// \param targs - target DAG nodes
  /// \param n_targs - number of target DAG nodes
  /// \param tint - internal target DAG nodes
  /// \param n_tint - number of internal target DAG nodes
  /// \param sint - internal source DAG nodes
  /// \param n_sint - number of internal source DAG nodes
  ///
  /// \returns - HPX_SUCCESS
  static int termination_detection_handler(hpx_addr_t done,
                                           std::vector<DAGNode *> *targs,
                                           size_t n_targs,
                                           std::vector<DAGNode *> *tint,
                                           size_t n_tint,
                                           std::vector<DAGNode *> *sint,
                                           size_t n_sint) {
    // If this is insufficiently parallel, we can always make this action
    // call itself with smaller and smaller chunks of the array.
    int nskip{0};

    int myrank = hpx_get_my_rank();
    for (size_t i = 0; i < n_targs; ++i) {
      assert((*targs)[i] != nullptr);
      if ((*targs)[i]->locality == myrank) {
        assert((*targs)[i]->global_addx != HPX_NULL);
        hpx_call_when((*targs)[i]->global_addx, done, hpx_lco_set_action,
                      HPX_NULL, nullptr, 0);
      } else {
        ++nskip;
      }
    }

    for (size_t i = 0; i < n_tint; ++i) {
      assert((*tint)[i] != nullptr);
      if ((*tint)[i]->locality == myrank) {
        assert((*tint)[i]->global_addx != HPX_NULL);
        hpx_call_when((*tint)[i]->global_addx, done, hpx_lco_set_action,
                      HPX_NULL, nullptr, 0);
      } else {
        ++nskip;
      }
    }

    for (size_t i = 0; i < n_sint; ++i) {
      assert((*sint)[i] != nullptr);
      if ((*sint)[i]->locality == myrank) {
        assert((*sint)[i]->global_addx != HPX_NULL);
        hpx_call_when((*sint)[i]->global_addx, done, hpx_lco_set_action,
                      HPX_NULL, nullptr, 0);
      } else {
        ++nskip;
      }
    }

    hpx_lco_and_set_num(done, nskip, HPX_NULL);

    return HPX_SUCCESS;
  }

  /// Action to destroy the DAG LCOs
  ///
  /// \param nodes - the DAG nodes
  /// \param n_nodes - the number of nodes
  ///
  /// \returns - HPX_SUCCESS
  static int destroy_DAG_LCOs_handler(DAGNode **nodes, size_t n_nodes,
                                      int type) {
    int myrank = hpx_get_my_rank();
    for (size_t i = 0; i < n_nodes; ++i) {
      if (nodes[i]->locality == myrank) {
        assert(nodes[i]->global_addx != HPX_NULL);
        if (type) {
          // NOTE: destroy() does not care about the second argument to the
          // constructor, so we put in some junk
          auto temp = targetlco_t{nodes[i]->global_addx, 0};
          temp.destroy();
        } else {
          auto temp = expansionlco_t{nodes[i]->global_addx};
          temp.destroy();
        }
      }
    }

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
  int *rank_map_;             /// map unif grid index to rank
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
  static hpx_action_t source_apply_method_;
  static hpx_action_t source_apply_method_child_done_;
  static hpx_action_t target_apply_method_;
  static hpx_action_t destroy_DAG_LCOs_;
  static hpx_action_t termination_detection_;
  static hpx_action_t create_S_expansions_from_DAG_;
  static hpx_action_t create_T_expansions_from_DAG_;
  static hpx_action_t edge_lists_;
  static hpx_action_t instigate_dag_eval_;
  static hpx_action_t instigate_dag_eval_remote_;
};

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::domain_geometry_init_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::domain_geometry_op_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::set_domain_geometry_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::init_partition_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::recv_points_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::send_points_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::create_dual_tree_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::finalize_partition_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::source_apply_method_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::source_apply_method_child_done_ =
    HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::target_apply_method_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::destroy_DAG_LCOs_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::termination_detection_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::create_S_expansions_from_DAG_ =
    HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::create_T_expansions_from_DAG_ =
    HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::edge_lists_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::instigate_dag_eval_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::instigate_dag_eval_remote_ =
    HPX_ACTION_NULL;


} // namespace dashmm


#endif
