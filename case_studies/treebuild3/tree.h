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


/// A node of the tree.
template <typename Record>
class Node {
 public:
  using record_t = Record;
  using node_t = Node<Record>;
  using arrayref_t = ArrayRef<Record>;

  Node() : idx{}, parts{}, parent{nullptr}, first_{0} {
    for (int i = 0; i < 8; ++i) {
      child[i] = nullptr;
    }
    sema_ = HPX_NULL;
    complete_ = HPX_NULL;
  }

  /// Constuct with a known index.
  Node(Index idx) : idx{idx}, parts{}, parent{nullptr}, first_{0} {
    for (int i = 0; i < 8; ++i) {
      child[i] = nullptr;
    }
    sema_ = HPX_NULL;
    complete_ = hpx_lco_and_new(8);
  }

  /// Construct with an index, a particle segment and a parent
  Node(Index idx, arrayref_t parts, node_t *parent)
      : idx{idx}, parts{parts}, parent{parent}, first_{0} {
    for (int i = 0; i < 8; ++i) {
      child[i] = nullptr;
    }
    sema_ = HPX_NULL;
    complete_ = hpx_lco_and_new(8);
  }

  ~Node() { }

  /// Returns the first record into which new records may be copied.
  size_t first() const {return first_;}

  /// Returns the LCO signaling completion of partitioning for this node
  hpx_addr_t complete() const {return complete_;}

  /// Gives the number of particles in the segment owned by this node
  size_t num_parts() const {return parts.n();}

  /// Create a completion detection and gate
  void add_completion() {
    assert(complete_ == HPX_NULL);
    complete_ = hpx_lco_and_new(8);
  }

  /// Create semaphore as needed
  void add_lock() {
    assert(sema_ == HPX_NULL);
    sema_ = hpx_lco_sema_new(1);
  }

  /// Delete semaphore
  void delete_lock() {
    assert(sema_ != HPX_NULL);
    hpx_lco_delete_sync(sema_);
    sema_ = HPX_NULL;
  }

  /// Lock the semaphore
  void lock() const {
    // TODO consider making this an if with graceful error case
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
  /// \param incr - the number of records that has just been added to the
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
  void partition(int threshold, DomainGeometry *geo) {
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
                   &cnd, &geo, &threshold);
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
    int count = 1;
    for (int i = 0; i < 8; ++i) {
      if (child[i] != nullptr) {
        count += child[i]->n_descendants();
      }
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
  /// \param branch - buffers holding the tree structure
  /// \param tree -
  /// \param parent - the index of the parent of this node
  /// \param curr - the current slot in the buffers; this is updated during
  ///               the call to compress()
  void compress(int *branch, int *tree, int parent, int &curr) const {
    for (int i = 0; i < 8; ++i) {
      if (child[i] != nullptr) {
        branch[curr] = i; // tracks which child exists
        tree[curr] = parent; // tracks the parent of the node being processed
        curr++; // Move onto the next slot
        // curr - 1 is the parent location for the subtree rooted at child[i]
        child[i]->compress(branch, tree, curr - 1, curr);
      }
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

  // This is not recursive because we are working with HPX-5's small stack
  // size probably.
  static void destroy_branch(node_t *root) {
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
  // TODO
  // DAGInfo dag;                    /// The DAG info for this node

  static hpx_action_t partition_node_;

 private:
  friend class NodeRegistrar<Record>;

  static int partition_node_handler(node_t *n, DomainGeometry *geo,
                                    int threshold) {
    n->partition(threshold, geo);
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
  using arrayref_t = ArrayRef<Record>;
  using tree_t = Tree<Source, Target, Record, Expansion, Method, DistroPolicy>;
  using dualtree_t = DualTree<Source, Target, Expansion, Method, DistroPolicy>;

  Tree() : root_{nullptr}, unif_grid_{nullptr}, unif_done_{HPX_NULL},
           sorted_{} { }

  // Some basic setup during initial tree construction
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

  // Tear down the data in the tree
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

  // This will assign the points to the uniform grid. This gives the points
  // the id (in the Morton Key sense) of the box to which they are assigned,
  // and it will count the numbers in each box.
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

  // This will rearrange the particles into their bin order. This is a stable
  // reordering.
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

  // This is the action on the other side that receives the partitioned tree.
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
        const int *branch = &compressed_tree[2];
        const int *tree = &compressed_tree[2 + n_nodes];
        grid[id].extract(branch, tree, n_nodes);
      }

      hpx_lco_and_set_num(grid[id].complete(), 8, HPX_NULL);
    } else {
      Node<Target> *grid = local_tree->unif_grid_target();

      if (n_nodes) {
        const int *branch = &compressed_tree[2];
        const int *tree = &compressed_tree[2 + n_nodes];
        grid[id].extract(branch, tree, n_nodes);
      }

      hpx_lco_and_set_num(grid[id].complete(), 8, HPX_NULL);
    }

    return HPX_SUCCESS;
  }

  // This action is called once the individual grids are done.
  // Also, this is where the points are finally rearranged. All other work has
  // been to set up their eventual index in the sorted situation. Here is the
  // actual shuffling.
  //
  // This will send one message for each grid box to all other localities.
  // This does allow for maximum parallelism. It is likely the right approach.
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
      int pos = 0;
      curr->compress(branch, tree, -1, pos);
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

  // This sets up some orgnaizational structures. The return value is an
  // array of offsets in the final set of points for each part of the uniform
  // grid. This is basically just setup. There is a chance that some work occurs
  // in one branch. Before leaving, this will bring the points that do not have
  // to change rank into their correct location. If it is detected that all of
  // the points for that part of the tree have arrived, then the partitioning
  // work will begin.
  static int *init_point_exchange(tree_t *tree, int first, int last,
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
          auto sorted = curr->parts.pin();
          memcpy(sorted.value() + curr->first(),
                 &temp[local_offset[i]],
                 sizeof(record_t) * local_count[i]);

          if (curr->increment_first(local_count[i])) {
            // This grid does not expect remote points.
            // Spawn adaptive partitioning
            hpx_call(HPX_HERE, node_t::partition_node_, HPX_NULL,
                     &curr, &geo, &threshold);
          }
        }
      }
    }

    // this is actually part of tree now, and not dual tree as is done here
    tree->sorted_ = sorted_ref;

    hpx_lco_and_set(tree->unif_done_, HPX_NULL);

    return global_offset;
  }

  // This action merges particular points with the sorted list. Also, if this
  // is the last set of points that are merged, this will go ahead and start
  // the adaptive partitioning of that part of the local tree.
  static int merge_points_handler(record_t *temp, node_t *n, int n_arrived,
                                  int n_total, hpx_addr_t rwgas) {
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
      hpx_call(HPX_HERE, node_t::partition_node_, HPX_NULL,
               &n, &geoarg, &thresh);
    }
    n->unlock();

    return HPX_SUCCESS;
  }

  static uint64_t split(unsigned k) {
    uint64_t split = k & 0x1fffff;
    split = (split | split << 32) & 0x1f00000000ffff;
    split = (split | split << 16) & 0x1f0000ff0000ff;
    split = (split | split << 8)  & 0x100f00f00f00f00f;
    split = (split | split << 4)  & 0x10c30c30c30c30c3;
    split = (split | split << 2)  & 0x1249249249249249;
    return split;
  }

  static uint64_t morton_key(unsigned x, unsigned y, unsigned z) {
    uint64_t key = 0;
    key |= split(x) | split(y) << 1 | split(z) << 2;
    return key;
  }


private:
  friend class DualTree<Source, Target, Expansion, Method, DistroPolicy>;
  friend class TreeRegistrar<Source, Target, Record, Expansion, Method,
                             DistroPolicy>;

  node_t *root_;
  node_t *unif_grid_;
  hpx_addr_t unif_done_;
  arrayref_t sorted_;

  static hpx_action_t setup_basics_;
  static hpx_action_t delete_tree_;
  static hpx_action_t recv_node_;
  static hpx_action_t send_node_;
  static hpx_action_t assign_points_;
  static hpx_action_t group_points_;
  static hpx_action_t merge_points_;
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


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


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

  DualTree()
    : domain_{}, refinement_limit_{1}, unif_level_{1}, dim3_{8},
      unif_count_{HPX_NULL}, unif_count_value_{nullptr},
      distribute_{nullptr}, method_{}, source_tree_{nullptr},
      target_tree_{nullptr} { }

  DualTree(const dualtree_t &other) = delete;
  dualtree_t &operator=(const dualtree_t &other) = delete;

  // TODO define move construction and assignment?
  // Is this needed?

  // simple accessors
  int *unif_count_src(size_t i = 0) const {return &unif_count_value_[i];}
  int *unif_count_tar(size_t i = 0) const {
    return &unif_count_value_[i + dim3_];
  }
  const DomainGeometry *domain() const {return &domain_;}
  int refinement_limit() const {return refinement_limit_;}

  sourcenode_t *unif_grid_source() {return source_tree_->unif_grid_;}

  targetnode_t *unif_grid_target() {return target_tree_->unif_grid_;}

  size_t sorted_src_count() const {return source_tree_->sorted_.n_tot();}
  size_t sorted_tar_count() const {return target_tree_->sorted_.n_tot();}
  sourcearraydata_t sorted_src() const {return source_tree_->sorted_.pin();}
  targetarraydata_t sorted_tar() const {return target_tree_->sorted_.pin();}

  const method_t &method() const {return method_;}
  void set_method(const method_t &method) {method_ = method;}

  // more complex things
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

  int first(int rank) const {return rank == 0 ? 0 : distribute_[rank - 1] + 1;}
  int last(int rank) const {return distribute_[rank];}


  // More DASHMM stuff
  // DAG create_DAG(bool same_sandt)
  // void collect_DAG_nodes(DAG &dag)
  // void create_expansions_from_DAG(int n_digits)
  // hpx_addr_t setup_termination_detection(DAG &dag)
  // void setup_edge_lists(DAG &dag)
  // void start_DAG_evaluation()
  // void destroy_DAG_LCOs(DAG &dag)


  // External interface - these are likely the most important

  // This should be called from inside HPX-5.
  //
  // Also, we want this to return as soon as the work is started. In this way,
  // we can do whatever overlap is possible. Then there should be some interface
  // to be sure it is complete or something. Possibly even another version of
  // create that is create_sync. This means whatever overlap stuff we have going
  // will have to be saved in the rankwise data.
  //
  // This is to be called from a single thread
  static RankWise<dualtree_t> create(int threshold, Array<Source> sources,
                                     Array<Target> targets) {
    hpx_addr_t domain_geometry = compute_domain_geometry(sources, targets);
    RankWise<dualtree_t> retval = setup_basic_data(threshold, domain_geometry);
    hpx_lco_delete_sync(domain_geometry);
    return retval;
  }

  // This should be called from inside HPX-5
  //
  // This is to be called from a single thread
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

  // This should be called from inside HPX-5
  //
  // This is to be called from a single thread
  static void destroy(RankWise<dualtree_t> global_tree) {
    hpx_addr_t rwtree = global_tree.data();
    hpx_bcast_rsync(finalize_partition_, &rwtree);

    auto tree = global_tree.here();
    hpx_lco_delete_sync(tree->unif_count_);
  }


 private:
  friend class Registrar<Source, Target, Expansion, Method, DistroPolicy>;

  static int set_domain_geometry_handler(hpx_addr_t sources_gas,
                                         hpx_addr_t targets_gas,
                                         hpx_addr_t domain_geometry) {
    Array<source_t> sources{sources_gas};
    sourceref_t src_ref = sources.ref();
    sourcearraydata_t src_data = src_ref.pin();
    Source *s = src_data.value();

    Array<target_t> targets{targets_gas};
    targetref_t trg_ref = targets.ref();
    targetarraydata_t trg_data = trg_ref.pin();
    Target *t = trg_data.value();

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

    for (size_t i = 0; i < trg_ref.n(); ++i) {
      var[0] = fmin(var[0], t[i].position.x());
      var[1] = fmax(var[1], t[i].position.x());
      var[2] = fmin(var[2], t[i].position.y());
      var[3] = fmax(var[3], t[i].position.y());
      var[4] = fmin(var[4], t[i].position.z());
      var[5] = fmax(var[5], t[i].position.z());
    }

    hpx_lco_set_lsync(domain_geometry, sizeof(double) * 6, var, HPX_NULL);

    return HPX_SUCCESS;
  }

  static void domain_geometry_init_handler(double *values,
                                           const size_t UNUSED) {
    values[0] = 1e50; // xmin
    values[1] = -1e50; // xmax
    values[2] = 1e50; // ymin
    values[3] = -1e50; // ymax
    values[4] = 1e50; // zmin
    values[5] = -1e50; // zmax
  }

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
                                            Array<Target> targets) {
    // Create a reduction LCO
    hpx_addr_t domain_geometry =
      hpx_lco_reduce_new(hpx_get_num_ranks(), sizeof(double) * 6,
                         domain_geometry_init_,
                         domain_geometry_op_);

    // Launch the reduction actions
    hpx_addr_t sglob = sources.data();
    hpx_addr_t tglob = targets.data();
    hpx_bcast_lsync(set_domain_geometry_, HPX_NULL,
                    &sglob, &tglob, &domain_geometry);

    return domain_geometry;
  }

  static int init_partition_handler(hpx_addr_t rwdata, hpx_addr_t count,
                                    int limit, hpx_addr_t domain_geometry) {
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

  static RankWise<dualtree_t> setup_basic_data(int threshold,
                                           hpx_addr_t domain_geometry) {
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
    hpx_bcast_rsync(init_partition_, &rwdata, &ucount, &threshold,
                    &domain_geometry);

    return retval;
  }

  // Given the global counts, this will partition the uniform grid among the
  // available localities. There is perhaps some room for simplification here.
  // I would rather shoot for fixed targets instead of aiming to take a fair
  // fraction of whatever is left. But there might not be any real performance
  // impact either way.
  //
  // TODO: Note that this would be a target for a second member of the
  // Distribution Policy
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

  // Get counting and sort the local points
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
    hpx_call(HPX_HERE, targettree_t::group_points_, group_done,
             &p_t, &n_targets, &tree->dim3_, &gid_of_targets, &local_tcount,
             &local_offset_t);
    hpx_lco_wait(group_done);
    hpx_lco_delete_sync(group_done);

    delete [] gid_of_sources;
    delete [] gid_of_targets;

    return local_count;
  }

  // This is the 'far-side' of the send points message. This action merges the
  // incoming points into the sorted list and will then spawn the adaptive
  // partition if this happens to be the last block for a given uniform grid.
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
        sourcenode_t *ns = &stree->unif_grid_[i];
        int incoming_ns = count_s[i - first];
        if (incoming_ns) {
          hpx_call(HPX_HERE, sourcetree_t::merge_points_, done,
                   &recv_s, &ns, &incoming_ns,
                   local_tree->unif_count_src(i), rwarg);
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
        targetnode_t *nt = &ttree->unif_grid_[i];
        int incoming_nt = count_t[i - first];
        if (incoming_nt) {
          hpx_call(HPX_HERE, targettree_t::merge_points_, done,
                   &recv_t, &nt, &incoming_nt,
                   local_tree->unif_count_tar(i), rwarg);
          recv_t += incoming_nt;
        } else {
          hpx_lco_and_set(done, HPX_NULL);
        }
      }
    } else {
      hpx_lco_and_set_num(done, range, HPX_NULL);
    }

    // Wait until the data has been merged before releasing the parcel
    // NOTE: This 'done' will trigger once points are merged. The action that
    // triggers this will spawn more work, but not in a synchronous way.
    hpx_lco_wait(done);
    hpx_lco_delete_sync(done);
    return HPX_SUCCESS;
  }

  // This is pretty straightforward. This rank will send to the given rank all
  // those particles that are to be shipped out. There is nothing too
  // complicated in this, just some indexing and so forth.
  //
  // NOTE: This would be a bit nicer looking using the Buffer types.
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
    int send_ns = 0, send_nt = 0;
    for (int i = first; i <= last; ++i) {
      send_ns += count_s[i];
      send_nt += count_t[i];
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

  // This appears to be the main action that creates the trees. This will
  // organize and call out to the other actions.
  // NOTE: One thing is sure, this ought to be factored
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
    int *global_offset_s = sourcetree_t::init_point_exchange(
        tree->source_tree_, firstarg, lastarg, tree->unif_count_src(),
        local_scount, local_offset_s, &tree->domain_, tree->refinement_limit_,
        p_s, ns);
    int *global_offset_t = targettree_t::init_point_exchange(
        tree->target_tree_, firstarg, lastarg, tree->unif_count_tar(),
        local_tcount, local_offset_t, &tree->domain_, tree->refinement_limit_,
        p_t, nt);

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
    hpx_addr_t old_tar_data = targets.replace(tree->target_tree_->sorted_);
    hpx_gas_free_sync(old_src_data);
    hpx_gas_free_sync(old_tar_data);

    delete [] local_count;
    delete [] local_offset_s;
    delete [] local_offset_t;
    delete [] global_offset_s;
    delete [] global_offset_t;

    return HPX_SUCCESS;
  }

  static int finalize_partition_handler(hpx_addr_t rwtree) {
    RankWise<dualtree_t> global_tree{rwtree};
    auto tree = global_tree.here();
    tree->clear_data();
    return HPX_SUCCESS;
  }

  //
  // Now for the data members
  //

  DomainGeometry domain_;
  int refinement_limit_;
  int unif_level_;
  int dim3_;
  hpx_addr_t unif_count_;
  int *unif_count_value_;
  int *distribute_;
  method_t method_;

  sourcetree_t *source_tree_;
  targettree_t *target_tree_;


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
