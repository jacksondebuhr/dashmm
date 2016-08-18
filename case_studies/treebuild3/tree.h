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
void registrar();
uint64_t morton_key(unsigned x, unsigned y, unsigned z);



/// A node of the tree.
class Node {
 public:
  Node() : idx{}, parts{}, parent{nullptr}, first_{0} {
    for (int i = 0; i < 8; ++i) {
      child[i] = nullptr;
    }
    sema_ = HPX_NULL;
    complete_ = HPX_NULL;
  }

  /// Constuct with a known index.
  Node(Index idx)
      : idx{idx}, parts{}, parent{nullptr}, first_{0} {
    for (int i = 0; i < 8; ++i) {
      child[i] = nullptr;
    }
    sema_ = HPX_NULL;
    complete_ = hpx_lco_and_new(8);
  }

  /// Construct with an index, a particle segment and a parent
  Node(Index idx, ArrayRef<Point> parts, Node *parent)
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

    if (parent) {
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
      Point *p = p_local.value();

      // Sort the particles among the children
      Point *splits[9]{};
      splits[0] = p;
      splits[8] = &p[num_points];

      auto z_comp = [&center_z](Point &a) {
        return a.z() < center_z;
      };
      splits[4] = std::partition(splits[0], splits[8], z_comp);

      auto y_comp = [&center_y](Point &a) {
        return a.y() < center_y;
      };
      splits[2] = std::partition(splits[0], splits[4], y_comp);
      splits[6] = std::partition(splits[4], splits[8], y_comp);

      auto x_comp = [&center_x](Point &a) {
        return a.x() < center_x;
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
          Node *cnd = new Node{idx.child(i), cparts, this};
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

    // TODO I think I don't want to do this one this way, favoring instead
    // a more typical structure. I will think about it
    Node *descendants = new Node[n_nodes];

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
      Node *curr = &descendants[i];
      Node *parent = (pos < 0 ? this : &descendants[pos]);

      curr->parent = parent;
      curr->idx = parent->idx.child(which);
      parent->child[which] = curr;
    }
  }

  /// Destroy the node - this will recursively destroy children, freeing up
  /// the two LCOs associated with the node. is allocated_in_array is true,
  /// this will not destroy the node itself, as that would cause trouble.
  void destroy(bool allocated_in_array) {
    for (int i = 0; i < 8; ++i) {
      if (child[i]) {
        child[i]->destroy(allocated_in_array);
      }
    }
    if (sema_ != HPX_NULL) {
      hpx_lco_delete_sync(sema_);
    }
    if (!allocated_in_array) {
      delete this;
    }
  }


  Index idx;                      /// index of the node
  ArrayRef<Point> parts;  /// segment for this node
  Node *parent;                           /// parent node
  Node *child[8];                         /// children of this node
  // TODO
  // DAGInfo dag;                            /// The DAG info for this node

  // TODO
  // We have to have this in node, and since stuff from DualTree needs to call
  // this we make it public. There could be some way to avoid the circular
  // reference needed to make this part of DualTree.
  static hpx_action_t partition_node_;

 private:
  friend void registrar();

  static int partition_node_handler(Node *n, DomainGeometry *geo,
                                    int threshold) {
    n->partition(threshold, geo);
    return HPX_SUCCESS;
  }

  size_t first_;            /// first record that is available
  hpx_addr_t sema_;         /// restrict concurrent modification
  hpx_addr_t complete_;     /// This is used to indicate that partitioning is
                            ///  complete
};


class DualTree {
 public:
  DualTree()
    : domain_{}, refinement_limit_{1}, unif_level_{1}, dim3_{8},
      unif_count_{HPX_NULL}, unif_count_value_{nullptr}, unif_grid_{nullptr},
      unif_done_{HPX_NULL}, distribute_{nullptr}, sorted_src_{},
      sorted_tar_{} { }

  DualTree(const DualTree &other) = delete;
  DualTree &operator=(const DualTree &other) = delete;

  // TODO define move construction and assignment?
  // Is this needed?

  // simple accessors and mutators
  // TODO this is a bit of a mess, is there a way to clean this up?
  //  When the various actions become part of the class, we can probably
  //  just hit these directly. This is way too much interface for the average
  //  user.
  int unif_level() const {return unif_level_;}
  int dim3() const {return dim3_;}
  Node *unif_grid_src(size_t i = 0) const {return &unif_grid_[i];}
  Node *unif_grid_tar(size_t i = 0) const {return &unif_grid_[dim3_ + i];}
  hpx_addr_t unif_count() const {return unif_count_;}
  int *unif_count_src(size_t i = 0) const {return &unif_count_value_[i];}
  int *unif_count_tar(size_t i = 0) const {
    return &unif_count_value_[i + dim3_];
  }
  hpx_addr_t unif_done() const {return unif_done_;}
  int refinement_limit() const {return refinement_limit_;}
  const DomainGeometry &domain() const {return domain_;}
  size_t sorted_src_count() const {return sorted_src_.n_tot();}
  size_t sorted_tar_count() const {return sorted_tar_.n_tot();}
  ArrayData<Point> sorted_src() const {return sorted_src_.pin();}
  ArrayData<Point> sorted_tar() const {return sorted_tar_.pin();}
  ArrayRef<Point> sorted_src_ref() const {return sorted_src_;}
  ArrayRef<Point> sorted_tar_ref() const {return sorted_tar_;}

  void set_unif_level(int l) {unif_level_ = l;}
  void set_dim3(int d) {dim3_ = d;}
  void set_unif_grid(Node *n) {unif_grid_ = n;}
  void set_unif_count(hpx_addr_t u) {unif_count_ = u;}
  void set_unif_count_value(int *u) {unif_count_value_ = u;}
  void set_unif_done(hpx_addr_t u) {unif_done_ = u;}
  void set_refinement_limit(int t) {refinement_limit_ = t;}
  void set_domain(const DomainGeometry &geo) {domain_ = geo;}
  void set_distribution(int *d) {
    assert(d != nullptr);
    distribute_ = d;
  }
  void set_sorted_src(ArrayRef<Point> s) {sorted_src_ = s;}
  void set_sorted_tar(ArrayRef<Point> t) {sorted_tar_ = t;}


  // Things to make it fit with DASHMM
  // const method_t &method() const {return method_;}
  // void set_method(const method_t &method) {method_ = method;}


  // more complex things
  void clear_data() {
    int rank = hpx_get_my_rank();

    delete [] unif_count_value_;
    hpx_lco_delete_sync(unif_done_);

    int b = first(rank);
    int e = last(rank);

    // NOTE: The difference here is that there are two different allocation
    // schemes for the nodes.

    // TODO: Add some parallelism here. Otherwise, this could take a long time.
    //  Most of the delay has been removed by making the lco delete themselves
    //  along the way, and only having those nodes that need semaphores have
    //  semaphores
    for (int i = 0; i < b; ++i) {
      Node *curr = &unif_grid_[i];
      curr->delete_lock();

      for (int j = 0; j < 8; ++j) {
        Node *child = curr->child[j];
        if (child) {
          child->destroy(true);
        }
      }

      for (int j = 0; j < 8; ++j) {
        Node *child = curr->child[j];
        if (child) {
          delete [] child;
          break;
        }
      }
    }

    for (int i = b; i <= e; ++i) {
      Node *curr = &unif_grid_[i];
      curr->delete_lock();

      for (int j = 0; j < 8; ++j) {
        Node *child = curr->child[j];
        if (child) {
          child->destroy(false);
        }
      }
    }

    for (int i = e + 1; i < dim3_; ++i) {
      Node *curr = &unif_grid_[i];
      curr->delete_lock();

      for (int j = 0; j < 8; ++j) {
        Node *child = curr->child[j];
        if (child) {
          child->destroy(true);
        }
      }

      for (int j = 0; j < 8; ++j) {
        Node *child = curr->child[j];
        if (child) {
          delete [] child;
          break;
        }
      }
    }

    for (int i = 0; i < b; ++i) {
      Node *curr = &unif_grid_[i + dim3_];
      curr->delete_lock();

      for (int j = 0; j < 8; ++j) {
        Node *child = curr->child[j];
        if (child) {
          child->destroy(true);
        }
      }

      for (int j = 0; j < 8; ++j) {
        Node *child = curr->child[j];
        if (child) {
          delete [] child;
          break;
        }
      }
    }

    for (int i = b; i <= e; ++i) {
      Node *curr = &unif_grid_[i + dim3_];
      curr->delete_lock();

      for (int j = 0; j < 8; ++j) {
        Node *child = curr->child[j];
        if (child) {
          child->destroy(false);
        }
      }
    }

    for (int i = e + 1; i < dim3_; ++i) {
      Node *curr = &unif_grid_[i + dim3_];
      curr->delete_lock();

      for (int j = 0; j < 8; ++j) {
        Node *child = curr->child[j];
        if (child) {
          child->destroy(true);
        }
      }

      for (int j = 0; j < 8; ++j) {
        Node *child = curr->child[j];
        if (child) {
          delete [] child;
          break;
        }
      }
    }

    delete [] unif_grid_;
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
  static RankWise<DualTree> create(int threshold, Array<Point> sources,
                                   Array<Point> targets) {
    hpx_addr_t domain_geometry = compute_domain_geometry(sources, targets);
    RankWise<DualTree> retval = setup_basic_data(threshold, domain_geometry);
    hpx_lco_delete_sync(domain_geometry);
    return retval;
  }

  // This should be called from inside HPX-5
  //
  // This is to be called from a single thread
  static hpx_addr_t partition(RankWise<DualTree> global_tree,
                              Array<Point> sources,
                              Array<Point> targets) {
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
  static void destroy(RankWise<DualTree> global_tree) {
    hpx_addr_t rwtree = global_tree.data();
    hpx_bcast_rsync(finalize_partition_, &rwtree);

    auto tree = global_tree.here();
    hpx_lco_delete_sync(tree->unif_count_);
  }


 private:
  friend void registrar();

  static int set_domain_geometry_handler(hpx_addr_t sources_gas,
                                         hpx_addr_t targets_gas,
                                         hpx_addr_t domain_geometry) {
    Array<Point> sources{sources_gas};
    ArrayRef<Point> src_ref = sources.ref();
    ArrayData<Point> src_data = src_ref.pin();
    Point *s = src_data.value();

    Array<Point> targets{targets_gas};
    ArrayRef<Point> trg_ref = targets.ref();
    ArrayData<Point> trg_data = trg_ref.pin();
    Point *t = trg_data.value();

    double var[6] = {1e50, -1e50, 1e50, -1e50, 1e50, -1e50};

    // TODO: add more parallelism
    for (size_t i = 0; i < src_ref.n(); ++i) {
      var[0] = fmin(var[0], s[i].x());
      var[1] = fmax(var[1], s[i].x());
      var[2] = fmin(var[2], s[i].y());
      var[3] = fmax(var[3], s[i].y());
      var[4] = fmin(var[4], s[i].z());
      var[5] = fmax(var[5], s[i].z());
    }

    for (size_t i = 0; i < trg_ref.n(); ++i) {
      var[0] = fmin(var[0], t[i].x());
      var[1] = fmax(var[1], t[i].x());
      var[2] = fmin(var[2], t[i].y());
      var[3] = fmax(var[3], t[i].y());
      var[4] = fmin(var[4], t[i].z());
      var[5] = fmax(var[5], t[i].z());
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
  static hpx_addr_t compute_domain_geometry(Array<Point> sources,
                                            Array<Point> targets) {
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
    RankWise<DualTree> global_tree{rwdata};
    auto tree = global_tree.here();

    int num_ranks = hpx_get_num_ranks();
    tree->unif_level_ = ceil(log(num_ranks) / log(8)) + 1;
    int dim = pow(2, tree->unif_level_);
    tree->dim3_ = pow(8, tree->unif_level_);
    tree->unif_count_ = count;
    tree->refinement_limit_ = limit;

    // We here allocate space for the result of the counting
    tree->unif_count_value_ = new int[tree->dim3_ * 2]();

    // Setup unif_done LCO
    tree->unif_done_ = hpx_lco_and_new(1);

    // Setup unif_grid
    Node *unif_grid = new Node[2 * tree->dim3_];
    for (int iz = 0; iz < dim; ++iz) {
      for (int iy = 0; iy < dim; ++iy) {
        for (int ix = 0; ix < dim; ++ix) {
          uint64_t mid = morton_key(ix, iy, iz);
          unif_grid[mid] = Node{Index{ix, iy, iz, tree->unif_level_}};
          unif_grid[mid + tree->dim3_] =
              Node{Index{ix, iy, iz, tree->unif_level_}};
          unif_grid[mid].add_lock();
          unif_grid[mid + tree->dim3_].add_lock();
        }
      }
    }
    tree->unif_grid_ = unif_grid;

    // Setup domain_
    double var[6];
    hpx_lco_get(domain_geometry, sizeof(double) * 6, &var);
    double length = fmax(var[1] - var[0],
                         fmax(var[3] - var[2], var[5] - var[4]));
    DomainGeometry geo{Point{(var[1] + var[0] - length) / 2,
                             (var[3] + var[2] - length) / 2,
                             (var[5] + var[4] - length) / 2}, length};
    tree->domain_ = geo;

    return HPX_SUCCESS;
  }

  static RankWise<DualTree> setup_basic_data(int threshold,
                                             hpx_addr_t domain_geometry) {
    RankWise<DualTree> retval{};
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

  // This will assign the points to the uniform grid. This gives the points
  // the id (in the Morton Key sense) of the box to which they are assigned,
  // and it will count the numbers in each box.
  static void assign_points_to_unif_grid(const Point *P, int npts,
                                         const DomainGeometry &geo,
                                         int unif_level, int *gid,
                                         int *count) {
    Point corner = geo.low();
    double scale = 1.0 / geo.size();

    // TODO: This is serial processing; is there some way to parallelize this?
    //   This would perhaps be worth timing.
    int dim = pow(2, unif_level);
    for (int i = 0; i < npts; ++i) {
      const Point *p = &P[i];
      int xid = std::min(dim - 1, (int)(dim * (p->x() - corner.x()) * scale));
      int yid = std::min(dim - 1, (int)(dim * (p->y() - corner.y()) * scale));
      int zid = std::min(dim - 1, (int)(dim * (p->z() - corner.z()) * scale));
      gid[i] = morton_key(xid, yid, zid);
      count[gid[i]]++;
    }
  }

  // This will rearrange the particles into their bin order. This is a stable
  // reordering.
  //
  // TODO: Note that this too is a serial operation. Could we do some on-rank
  // parallelism here?
  //
  // TODO: Work out what the correct way to parallelize this would be. And then
  // make that happen.
  static int *group_points_on_unif_grid(Point *p_in, int npts, int dim3,
                                        int *gid_of_points, const int *count) {
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
    Point source, save;
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

    return offset;
  }

  // New factor
  static int *sort_local_points(DualTree *tree, Point *p_s,
                                int n_sources, Point *p_t, int n_targets,
                                int **local_offset_s, int **local_offset_t) {
    int *local_count = new int[tree->dim3_ * 2]();
    int *local_scount = local_count;
    int *local_tcount = &local_count[tree->dim3_];

    int *gid_of_sources = new int[n_sources]();
    int *gid_of_targets = new int[n_targets]();
    // TODO: perhaps these are actions? That is a coarse parallelism - then
    // perhaps more might be added inside these functions
    assign_points_to_unif_grid(p_s, n_sources, tree->domain(),
                               tree->unif_level_, gid_of_sources,
                               local_scount);
    assign_points_to_unif_grid(p_t, n_targets, tree->domain(),
                               tree->unif_level_, gid_of_targets,
                               local_tcount);

    // Exchange counts
    hpx_lco_set(tree->unif_count_, sizeof(int) * tree->dim3_ * 2, local_count,
                HPX_NULL, HPX_NULL);

    // Put points of the same grid together while waiting for
    // counting to complete
    // TODO: Perhaps start these as actions to get at least that coarse
    // parallelism.
    *local_offset_s = group_points_on_unif_grid(p_s, n_sources, tree->dim3_,
                                                gid_of_sources, local_scount);
    *local_offset_t = group_points_on_unif_grid(p_t, n_targets, tree->dim3_,
                                                gid_of_targets, local_tcount);
    delete [] gid_of_sources;
    delete [] gid_of_targets;

    return local_count;
  }

  // This sets up some orgnaizational structures. The return value is an
  // array of offsets in the final set of points for each part of the uniform
  // grid. This is basically just setup. There is a chance that some work occurs
  // in one branch. Before leaving, this will bring the points that do not have
  // to change rank into their correct location. If it is detected that all of
  // the points for that part of the tree have arrived, then the partitioning
  // work will begin.
  static int *init_point_exchange(int rank, DualTree *tree,
                                  const int *local_count,
                                  const int *local_offset,
                                  const Point *temp, Node *n,
                                  char type) {
    int first = tree->first(rank);
    int last = tree->last(rank);
    int range = last - first + 1;

    ArrayRef<Point> sorted_ref{};
    int *global_count{nullptr};

    if (type == 's') {
      global_count = tree->unif_count_src();
    } else {
      global_count = tree->unif_count_tar();
    }

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
          1, sizeof(Point) * num_points, 0);
      assert(sorted_gas != HPX_NULL);
      sorted_ref = ArrayRef<Point>{sorted_gas, num_points, num_points};

      for (int i = first; i <= last; ++i) {
        Node *curr = &n[i];
        curr->parts = sorted_ref.slice(global_offset[i - first],
                                       global_count[i]);

        if (local_count[i]) {
          // Copy local points before merging remote points
          auto sorted = curr->parts.pin();
          memcpy(sorted.value() + curr->first(),
                 &temp[local_offset[i]],
                 sizeof(Point) * local_count[i]);

          // probably remove the else branch and just put the increment inside
          // the conditional here
          if (curr->increment_first(local_count[i])) {
            // This grid does not expect remote points.
            // Spawn adaptive partitioning
            // TODO: This is a mild cheat for the type...
            const DomainGeometry *arg = &(tree->domain_);
            int threshold = tree->refinement_limit_;
            hpx_call(HPX_HERE, Node::partition_node_, HPX_NULL,
                     &curr, &arg, &threshold);
          }
        }
      }
    }

    if (type == 's') {
      tree->sorted_src_ = sorted_ref;
    } else {
      tree->sorted_tar_ = sorted_ref;
    }

    return global_offset;
  }

  // This action merges particular points with the sorted list. Also, if this
  // is the last set of points that are merged, this will go ahead and start
  // the adaptive partitioning of that part of the local tree.
  static int merge_points_handler(Point *temp, Node *n,
                           int n_arrived, int n_total, char type,
                           hpx_addr_t rwgas) {
    RankWise<DualTree> global_tree{rwgas};
    auto local_tree = global_tree.here();

    // Note: all the pointers are local to the calling rank.
    n->lock();
    size_t first = n->first();
    auto localp = n->parts.pin();
    Point *p = localp.value();
    memcpy(p + first, temp, sizeof(Point) * n_arrived);

    if (n->increment_first(n_arrived)) {
      const DomainGeometry *geoarg = &(local_tree->domain_);
      int thresh = local_tree->refinement_limit_;

      if (type == 's') {
        hpx_call(HPX_HERE, Node::partition_node_, HPX_NULL,
                 &n, &geoarg, &thresh);
      } else {
        hpx_call(HPX_HERE, Node::partition_node_, HPX_NULL,
                 &n, &geoarg, &thresh);
      }
    }
    n->unlock();

    return HPX_SUCCESS;
  }

  // This is the 'far-side' of the send points message. This action merges the
  // incoming points into the sorted list and will then spawn the adaptive
  // partition if this happens to be the last block for a given uniform grid.
  static int recv_points_handler(void *args, size_t UNUSED) {
    hpx_addr_t *rwarg = static_cast<hpx_addr_t *>(args);
    RankWise<DualTree> global_tree{*rwarg};
    auto local_tree = global_tree.here();

    // Wait until the buffer is allocated before merging incoming messages
    // We could do this as a call when, but then we need to be aware of the
    // addresses for every rank's LCO. For now, we do this, as it is simpler.
    hpx_lco_wait(local_tree->unif_done_);

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
    Point *recv_s =
      reinterpret_cast<Point *>(static_cast<char *>(args) + sizeof(int) * 4 +
                                sizeof(hpx_addr_t) +
                                sizeof(int) * range * (recv_ns > 0) +
                                sizeof(int) * range * (recv_nt > 0));
    Point *recv_t = recv_s + recv_ns;

    hpx_addr_t done = hpx_lco_and_new(range * 2);

    if (recv_ns) {
      char type = 's';
      for (int i = first; i <= last; ++i) {
        Node *ns = local_tree->unif_grid_src(i);
        int incoming_ns = count_s[i - first];
        if (incoming_ns) {
          hpx_call(HPX_HERE, merge_points_, done,
                   &recv_s, &ns, &incoming_ns,
                   local_tree->unif_count_src(i), &type, rwarg);
          recv_s += incoming_ns;
        } else {
          hpx_lco_and_set(done, HPX_NULL);
        }
      }
    } else {
      hpx_lco_and_set_num(done, range, HPX_NULL);
    }

    if (recv_nt) {
      char type = 't';
      for (int i = first; i <= last; ++i) {
        Node *nt = local_tree->unif_grid_tar(i);
        int incoming_nt = count_t[i - first];
        if (incoming_nt) {
          hpx_call(HPX_HERE, merge_points_, done,
                   &recv_t, &nt, &incoming_nt,
                   local_tree->unif_count_tar(i),
                   &type, rwarg);
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
  // those particles that are to be shipped out. There is nothing too complicated
  // in this, just some indexing and so forth.
  //
  // NOTE: This would be a bit nicer looking using the Buffer types.
  static int send_points_handler(int rank, int *count_s, int *count_t,
                                 int *offset_s, int *offset_t,
                                 Point *sources, Point *targets,
                                 hpx_addr_t rwaddr) {
    RankWise<DualTree> global_tree{rwaddr};
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
      bytes += sizeof(int) * range + sizeof(Point) * send_ns;
    }
    if (send_nt) {
      bytes += sizeof(int) * range + sizeof(Point) * send_nt;
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
      memcpy(meta_s, &sources[offset_s[first]], sizeof(Point) * send_ns);
    }

    char *meta_t = meta_s + send_ns * sizeof(Point);
    if (send_nt) {
      memcpy(meta_t, &targets[offset_t[first]], sizeof(Point) * send_nt);
    }

    hpx_parcel_set_target(p, HPX_THERE(rank));
    hpx_parcel_set_action(p, recv_points_);
    hpx_parcel_send(p, HPX_NULL);

    return HPX_SUCCESS;
  }

  // This is the action on the other side that receives the partitioned tree.
  static int recv_node_handler(char *message_buffer, size_t UNUSED) {
    hpx_addr_t *rwdata = reinterpret_cast<hpx_addr_t *>(message_buffer);
    int *compressed_tree = reinterpret_cast<int *>(
                                message_buffer + sizeof(hpx_addr_t));
    RankWise<DualTree> global_tree{*rwdata};
    auto local_tree = global_tree.here();
    int type = compressed_tree[0];
    int id = compressed_tree[1];
    int n_nodes = compressed_tree[2];
    Node *curr{nullptr};
    if (type) {
      curr = local_tree->unif_grid_src(id);
    } else {
      curr = local_tree->unif_grid_tar(id);
    }

    if (n_nodes) {
      const int *branch = &compressed_tree[3];
      const int *tree = &compressed_tree[3 + n_nodes];
      curr->extract(branch, tree, n_nodes);
    }

    hpx_lco_and_set_num(curr->complete(), 8, HPX_NULL);

    return HPX_SUCCESS;
  }

  // This action is called once the individual grids are done.
  // Also, this is where the points are finally rearranged. All other work has
  // been to set up their eventual index in the sorted situation. Here is the
  // actual shuffling.
  //
  // This will send one message for each grid box to all other localities.
  // This does allow for maximum parallelism. It is likely the right approach.
  static int send_node_handler(Node *n, Point *sorted, int id, int type,
                               hpx_addr_t rwaddr) {
    Node *curr = &n[id];

    // TODO: I think this is okay to remove; check this
    //RankWise<DualTree> global_tree{rwaddr};
    //auto local_tree = global_tree.here();

    // Exclude curr as it is already allocated on remote localities
    int n_nodes = curr->n_descendants() - 1;
    size_t msgsize = sizeof(int) * (3 + n_nodes * 2) + sizeof(hpx_addr_t);
    char *message_buffer = new char[msgsize];
    hpx_addr_t *rwdata = reinterpret_cast<hpx_addr_t *>(message_buffer);
    int *compressed_tree = reinterpret_cast<int *>(
                                message_buffer + sizeof(hpx_addr_t));

    *rwdata = rwaddr;
    compressed_tree[0] = type; // source tree is 0, target tree is 1
    compressed_tree[1] = id; // where to merge
    compressed_tree[2] = n_nodes; // # of nodes
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

  // This appears to be the main action that creates the trees. This will
  // organize and call out to the other actions.
  // NOTE: One thing is sure, this ought to be factored
  static int create_dual_tree_handler(hpx_addr_t rwtree,
                                      hpx_addr_t sources_gas,
                                      hpx_addr_t targets_gas) {
    int rank = hpx_get_my_rank();
    int num_ranks = hpx_get_num_ranks();

    RankWise<DualTree> global_tree{rwtree};
    auto tree = global_tree.here();

    Array<Point> sources{sources_gas};
    ArrayRef<Point> src_ref = sources.ref(rank);
    ArrayData<Point> src_data = src_ref.pin();
    Point *p_s = src_data.value();
    int n_sources = src_ref.n();

    Array<Point> targets{targets_gas};
    ArrayRef<Point> trg_ref = targets.ref(rank);
    ArrayData<Point> trg_data = trg_ref.pin();
    Point *p_t = trg_data.value();
    int n_targets = trg_ref.n();

    // Assign points to uniform grid
    int *local_offset_s{nullptr};
    int *local_offset_t{nullptr};
    int *local_count = DualTree::sort_local_points(&*tree, p_s, n_sources,
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
    Node *ns = tree->unif_grid_src();
    Node *nt = tree->unif_grid_tar();

    int *global_offset_s = init_point_exchange(rank, &*tree, local_scount,
                                               local_offset_s, p_s, ns, 's');
    int *global_offset_t = init_point_exchange(rank, &*tree, local_tcount,
                                               local_offset_t, p_t, nt, 't');
    hpx_lco_and_set(tree->unif_done_, HPX_NULL);

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

    for (int r = 0; r < num_ranks; ++r) {
      int first = tree->first(r);
      int last = tree->last(r);
      int s{0}, t{1};

      if (r == rank) {
        for (int i = first; i <= last; ++i) {
          if (*(tree->unif_count_src(i)) == 0) {
            hpx_lco_and_set(dual_tree_complete, HPX_NULL);
          } else {
            Point *arg = tree->sorted_src().value();
            hpx_call_when_with_continuation(ns[i].complete(), HPX_HERE,
                                            send_node_, dual_tree_complete,
                                            hpx_lco_set_action, &ns, &arg,
                                            &i, &s, &rwtree);
            hpx_call_when(ns[i].complete(),
                          ns[i].complete(), hpx_lco_delete_action,
                          HPX_NULL, nullptr, 0);
          }

          if (*(tree->unif_count_tar(i)) == 0) {
            hpx_lco_and_set(dual_tree_complete, HPX_NULL);
          } else {
            Point *arg = tree->sorted_tar().value();
            hpx_call_when_with_continuation(nt[i].complete(), HPX_HERE,
                                            send_node_, dual_tree_complete,
                                            hpx_lco_set_action, &nt, &arg,
                                            &i, &t, &rwtree);
            hpx_call_when(nt[i].complete(),
                          nt[i].complete(), hpx_lco_delete_action,
                          HPX_NULL, nullptr, 0);
          }
        }
      } else {
        for (int i = first; i <= last; ++i) {
          if (*(tree->unif_count_src(i)) == 0) {
            hpx_lco_and_set(dual_tree_complete, HPX_NULL);
          } else {
            hpx_call_when_with_continuation(ns[i].complete(),
                dual_tree_complete, hpx_lco_set_action,
                ns[i].complete(), hpx_lco_delete_action,
                nullptr, 0);
          }

          if (*(tree->unif_count_tar(i)) == 0) {
            hpx_lco_and_set(dual_tree_complete, HPX_NULL);
          } else {
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
    hpx_addr_t old_src_data = sources.replace(tree->sorted_src_ref());
    hpx_addr_t old_tar_data = targets.replace(tree->sorted_tar_ref());
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
    RankWise<DualTree> global_tree{rwtree};
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
  Node *unif_grid_;
  hpx_addr_t unif_done_;

  int *distribute_;

  ArrayRef<Point> sorted_src_;
  ArrayRef<Point> sorted_tar_;

  static hpx_action_t domain_geometry_init_;
  static hpx_action_t domain_geometry_op_;
  static hpx_action_t set_domain_geometry_;
  static hpx_action_t init_partition_;
  static hpx_action_t merge_points_;
  static hpx_action_t recv_points_;
  static hpx_action_t send_points_;
  static hpx_action_t recv_node_;
  static hpx_action_t send_node_;
  static hpx_action_t create_dual_tree_;
  static hpx_action_t finalize_partition_;
};


#endif
