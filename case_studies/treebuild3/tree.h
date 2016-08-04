#ifndef __TREE_H__
#define __TREE_H__

#include <vector>
#include "hpx/hpx.h"

/// The main action
extern hpx_action_t main_action;

/// Meta data for the array.
///
/// The array is a cyclic allocation of ArrayMetaData objects, one per locality
/// that inside point to the local data (data). This gives the size and the
/// count of the records as well as the pointer to local memory.
struct ArrayMetaData {
  size_t size;
  int count;
  char *data;
};

/// Serves as index to indentify nodes without having to worry over floating
/// point artihmetic.
class Index {
 public:
  Index(): level_{0}, idx_{0, 0, 0} {}
  Index(int level, int x, int y, int z) : level_{level}, idx_{x, y, z} {}

  int level() const {return level_;}
  int x() const {return idx_[0];}
  int y() const {return idx_[1];}
  int z() const {return idx_[2];}

  Index parent(int num = 1) const {
    return Index{level_ - 1, idx_[0] >> num, idx_[1] >> num, idx_[2] >> num};
  }

  Index child(int which) const {
    return Index{level_ + 1, (idx_[0] << 1) + ((170 >> which) & 1),
        (idx_[1] << 1) + ((204 >> which) & 1),
        (idx_[2] << 1) + ((240 >> which) & 1)};
  }

 private:
  int level_;
  int idx_[3];
};

/// Convenience point class to represent a location
class Point {
 public:
  Point(): pos_{0.0, 0.0, 0.0} {}
  Point(double x, double y, double z): pos_{x, y, z} {}
  Point(const Point &pt) : pos_{pt.x(), pt.y(), pt.z()} {}

  double x() const {return pos_[0];}
  double y() const {return pos_[1];}
  double z() const {return pos_[2];}

 private:
  double pos_[3];
};

/// A node of the tree.
class Node {
 public:
  /// Default construct. Nothing much here
  Node() : idx_{}, first_{-1}, last_{-1}, parent_{nullptr} {
    for (int i = 0; i < 8; ++i)
      child_[i] = nullptr;
    sema_ = HPX_NULL;
    complete_ = HPX_NULL;
  }

  /// Constuct with a known index.
  Node(Index idx): idx_{idx}, first_{-1}, last_{-1}, parent_{nullptr} {
    for (int i = 0; i < 8; ++i)
      child_[i] = nullptr;
    sema_ = hpx_lco_sema_new(1);
    complete_ = hpx_lco_and_new(8);
  }

  Node(Index idx, int first, int last, Node *parent)
      : idx_{idx}, first_{first}, last_{last}, parent_{parent} {
    for (int i = 0; i < 8; ++i)
      child_[i] = nullptr;
    sema_ = hpx_lco_sema_new(1);
    complete_ = hpx_lco_and_new(8);
  }

  ~Node() { }

  /// These are all simple accessors
  Index index() const {return idx_;}
  int first() const {return first_;}
  int last() const {return last_;}
  Node *parent() const {return parent_;}
  Node *child(int i) const {return child_[i];}
  hpx_addr_t sema() const {return sema_;}
  hpx_addr_t complete() const {return complete_;}

  // These are simple mutators
  void set_index(Index idx) {idx_ = idx;}
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
  void partition(Point *P, int *swap, int *bin, int *map, int threshold,
                 double corner_x, double corner_y, double corner_z,
                 double size);

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
  Index idx_;               /// index of the node
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

