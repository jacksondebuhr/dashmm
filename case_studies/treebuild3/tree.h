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
  /// Default construct. Nothing much here
  Node() : idx_{}, first_{-1}, last_{-1}, parent_{nullptr} {
    for (int i = 0; i < 8; ++i) {
      child_[i] = nullptr;
    }
    sema_ = HPX_NULL;
    complete_ = HPX_NULL;
  }

  /// Constuct with a known index.
  Node(dashmm::Index idx)
      : idx_{idx}, first_{-1}, last_{-1}, parent_{nullptr} {
    for (int i = 0; i < 8; ++i) {
      child_[i] = nullptr;
    }
    sema_ = hpx_lco_sema_new(1);
    complete_ = hpx_lco_and_new(8);
  }

  Node(dashmm::Index idx, int first, int last, Node *parent)
      : idx_{idx}, first_{first}, last_{last}, parent_{parent} {
    for (int i = 0; i < 8; ++i) {
      child_[i] = nullptr;
    }
    sema_ = hpx_lco_sema_new(1);
    complete_ = hpx_lco_and_new(8);
  }

  ~Node() { }

  /// These are all simple accessors
  dashmm::Index index() const {return idx_;}
  int first() const {return first_;}
  int last() const {return last_;}
  Node *parent() const {return parent_;}
  Node *child(int i) const {return child_[i];}
  hpx_addr_t sema() const {return sema_;}
  hpx_addr_t complete() const {return complete_;}

  // These are simple mutators
  void set_index(dashmm::Index idx) {idx_ = idx;}
  void set_first(int first) {first_ = first;}
  void set_last(int last) {last_ = last;}
  void set_parent(Node *parent) {parent_ = parent;}
  void set_child(int i, Node *child) {child_[i] = child;}
  void set_sema(hpx_addr_t sema) {sema_ = sema;}
  void set_complete(hpx_addr_t complete) {complete_ = complete;}

  /// Partition the given points
  ///
  /// P - the point data; this is a per locality array in local memory
  /// swap - some temporary storage used in sorting the point data
  /// bin - temporary storage used in sorting the point data
  /// map -  a map indicating which index a given point will be in the sorted
  ///        array of points
  /// threshold - the partitioning threshold
  /// corner_x, corner_y, corner_z - position of the corner of this(?) node
  /// size - the size of this node
  ///
  /// TODO: we can likely merge the bin and swap arrays to save on memory
  void partition(dashmm::Point *P, int *swap, int *bin, int *map,
                 int threshold, dashmm::Point corner, double size);

  /// Return the number of descendants of this node - this is recursive, and
  /// collects all of them
  int n_descendants() const;


  void compress(int *branch, int *tree, int parent, int &curr) const;
  void extract(const int *branch, const int *tree, int n_nodes);

  /// Destroy the node - this will recursively destroy children, freeing up
  /// the two LCOs associated with the node. is allocated_in_array is true,
  /// this will not destroy the node itself, as that would cause trouble.
  void destroy(bool allocated_in_array);

 private:
  dashmm::Index idx_;       /// index of the node
  int first_;               /// first index into the local share of points that
                            ///  are represented by this node
  int last_;                /// last index into the local share of points that
                            ///  are represented by this node
  Node *parent_;            /// parent node
  Node *child_[8];          /// children of this node
  hpx_addr_t sema_;         /// UNKNOWN restrict concurrent modification?
  hpx_addr_t complete_;     /// This is used to indicate that partitioning is
                            ///  complete
};


class DualTree {
 public:
  DualTree()
    : domain_{}, threshold_{1}, unif_level_{1}, dim_3{8},
      unif_count_{HPX_NULL}, unif_count_value_{nullptr}, unif_grid_{nullptr},
      unif_done_{HPX_NULL}, distribute_{nullptr}, sorted_src_count_{0},
      sorted_src_{nullptr}, sorted_tar_count_{0}, sorted_tar_{nullptr},
      swap_src_{nullptr}, bin_src_{nullptr}, map_src_{nullptr},
      swap_tar_{nullptr}, bin_tar_{nullptr}, map_tar_{nullptr} { }

  // simple accessors and mutators
  int unif_level() const {return unif_level_;}
  int dim3() const {return dim3_;}
  hpx_addr_t unif_count() const {return unif_count_;}
  int threshold() const {return threshold_;}
  int *unif_count_value() const {return unif_count_value_;}
  hpx_addr_t unif_done() const {return unif_done_;}
  Node *unif_grid() const {return unif_grid_;}
  const DomainGeometry &domain() const {return domain_;}

  void set_unif_level(int l) {unif_level_ = l;}
  void set_dim3(int d) {dim3_ = d;}
  void set_unif_count(hpx_addr_t u) {unif_count_ = u;}
  void set_threshold(int t) {threshold_ = t;}
  void set_unif_count_value(int *u) {unif_count_value_ = u;}
  void set_unif_done(hpx_addr_t u) {unif_done_ = u;}
  void set_unif_grid(Node *n) {unif_grid_ = n;}
  void set_domain(const DomainGeometry &geo) {domain_ = geo;}

  // more complex things
  void clear_data();

 private:
  DomainGeometry domain_;
  int threshold_;

  int unif_level_;
  int dim3_;
  hpx_addr_t unif_count_;
  int *unif_count_value_;
  Node *unif_grid_;
  hpx_addr_t unif_done_;

  int *distribute_;

  size_t sorted_src_count_;
  Point *sorted_src_;
  size_t sorted_tar_count_;
  Point *sorted_tar_;

  int *swap_src_;
  int *bin_src_;
  int *map_src_;
  int *swap_tar_;
  int *bin_tar_;
  int *map_tar_;
}


// Interface that has been introduced
void create_dual_tree_old(dashmm::Array<dashmm::Point> &sources,
                      dashmm::Array<dashmm::Point> &targets);


// TODO: Are these eventually put into DualTree itself?
RankWise<DualTree> dual_tree_create(int threshold,
                                    dashmm::Array<dashmm::Point> sources,
                                    dashmm::Array<dashmm::Point> targets);
void dual_tree_destroy(RankWise<DualTree> global_tree);

#endif
