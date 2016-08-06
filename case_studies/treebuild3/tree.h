#ifndef __TREE_H__
#define __TREE_H__

#include <vector>
#include "hpx/hpx.h"

#include "dashmm/array.h"
#include "dashmm/index.h"
#include "dashmm/point.h"


// Interface that has been introduced
void compute_domain_geometry(dashmm::Array<dashmm::Point> sources,
                             dashmm::Array<dashmm::Point> targets);
void setup_basic_data(int threshold);
void create_dual_tree(dashmm::Array<dashmm::Point> &sources,
                      dashmm::Array<dashmm::Point> &targets);
void finalize_partition();


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

#endif

