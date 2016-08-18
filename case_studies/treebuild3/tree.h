#ifndef __TREE_H__
#define __TREE_H__

#include <vector>
#include "hpx/hpx.h"

#include "dashmm/array.h"
#include "dashmm/index.h"
#include "dashmm/point.h"

#include "rankwise.h"


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
  Node(dashmm::Index idx)
      : idx{idx}, parts{}, parent{nullptr}, first_{0} {
    for (int i = 0; i < 8; ++i) {
      child[i] = nullptr;
    }
    sema_ = HPX_NULL;
    complete_ = hpx_lco_and_new(8);
  }

  /// Construct with an index, a particle segment and a parent
  Node(dashmm::Index idx, dashmm::ArrayRef<dashmm::Point> parts, Node *parent)
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

  void delete_lock() {
    assert(sema_ != HPX_NULL);
    hpx_lco_delete_sync(sema_);
    sema_ = HPX_NULL;
  }

  void lock() const {
    // TODO consider making this an if with graceful error case
    assert(sema_ != HPX_NULL);
    hpx_lco_sema_p(sema_);
  }

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
  void partition(int threshold, dashmm::DomainGeometry *geo);

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


  dashmm::Index idx;                      /// index of the node
  dashmm::ArrayRef<dashmm::Point> parts;  /// segment for this node
  Node *parent;                           /// parent node
  Node *child[8];                         /// children of this node

 private:
  size_t first_;            /// first record that is available
  hpx_addr_t sema_;         /// restrict concurrent modification
  hpx_addr_t complete_;     /// This is used to indicate that partitioning is
                            ///  complete
};


class DualTree {
 public:
  DualTree()
    : domain_{}, threshold_{1}, unif_level_{1}, dim3_{8},
      unif_count_{HPX_NULL}, unif_count_value_{nullptr}, unif_grid_{nullptr},
      unif_done_{HPX_NULL}, distribute_{nullptr}, sorted_src_{},
      sorted_tar_{} { }

  // simple accessors and mutators
  int unif_level() const {return unif_level_;}
  int dim3() const {return dim3_;}
  hpx_addr_t unif_count() const {return unif_count_;}
  int threshold() const {return threshold_;}
  int *unif_count_value() const {return unif_count_value_;}
  hpx_addr_t unif_done() const {return unif_done_;}
  Node *unif_grid() const {return unif_grid_;}
  const dashmm::DomainGeometry &domain() const {return domain_;}
  int *distribute() const {return distribute_;}
  size_t sorted_src_count() const {return sorted_src_.n_tot();}
  size_t sorted_tar_count() const {return sorted_tar_.n_tot();}
  dashmm::ArrayData<dashmm::Point> sorted_src() const {
    return sorted_src_.pin();
  }
  dashmm::ArrayData<dashmm::Point> sorted_tar() const {
    return sorted_tar_.pin();
  }
  dashmm::ArrayRef<dashmm::Point> sorted_src_ref() const {
    return sorted_src_;
  }
  dashmm::ArrayRef<dashmm::Point> sorted_tar_ref() const {
    return sorted_tar_;
  }

  void set_unif_level(int l) {unif_level_ = l;}
  void set_dim3(int d) {dim3_ = d;}
  void set_unif_count(hpx_addr_t u) {unif_count_ = u;}
  void set_threshold(int t) {threshold_ = t;}
  void set_unif_count_value(int *u) {unif_count_value_ = u;}
  void set_unif_done(hpx_addr_t u) {unif_done_ = u;}
  void set_unif_grid(Node *n) {unif_grid_ = n;}
  void set_domain(const dashmm::DomainGeometry &geo) {domain_ = geo;}
  void set_distribute(int *d) {distribute_ = d;}
  void set_sorted_src(dashmm::ArrayRef<dashmm::Point> s) {sorted_src_ = s;}
  void set_sorted_tar(dashmm::ArrayRef<dashmm::Point> t) {sorted_tar_ = t;}

  // more complex things
  void clear_data();
  int first(int rank) const {return rank == 0 ? 0 : distribute_[rank - 1] + 1;}
  int last(int rank) const {return distribute_[rank];}

 private:
  dashmm::DomainGeometry domain_;
  int threshold_;

  int unif_level_;
  int dim3_;
  hpx_addr_t unif_count_;
  int *unif_count_value_;
  Node *unif_grid_;
  hpx_addr_t unif_done_;

  int *distribute_;

  dashmm::ArrayRef<dashmm::Point> sorted_src_;
  dashmm::ArrayRef<dashmm::Point> sorted_tar_;
};


// TODO: Are these eventually put into DualTree itself?
RankWise<DualTree> dual_tree_create(int threshold,
                                    dashmm::Array<dashmm::Point> sources,
                                    dashmm::Array<dashmm::Point> targets);
hpx_addr_t dual_tree_partition(RankWise<DualTree> global_tree,
                               dashmm::Array<dashmm::Point> sources,
                               dashmm::Array<dashmm::Point> targets);
void dual_tree_destroy(RankWise<DualTree> global_tree);


#endif
