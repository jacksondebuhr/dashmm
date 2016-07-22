#ifndef __TREE_H__
#define __TREE_H__

#include <vector>
#include "hpx/hpx.h"

extern hpx_action_t main_action; 

struct ArrayMetaData {
  size_t size; 
  int count; 
  char *data;
}; 

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

class Node {
public:
  Node() : idx_{}, first_{0}, last_{0}, parent_{nullptr} {
    for (int i = 0; i < 8; ++i) 
      child_[i] = nullptr; 
  }

  Node(Index idx): idx_{idx}, first_{0}, last_{0}, parent_{nullptr} {
    for (int i = 0; i < 8; ++i) 
      child_[i] = nullptr; 
    sema_ = hpx_lco_sema_new(1); 
    complete_ = hpx_lco_and_new(8); 
  }

  Node(Index idx, int first, int last, Node *parent) : 
    idx_{idx}, first_{first}, last_{last}, parent_{parent}
  {
    for (int i = 0; i < 8; ++i) 
      child_[i] = nullptr; 
    sema_ = hpx_lco_sema_new(1); 
    complete_ = hpx_lco_and_new(8); 
  }

  ~Node() {}

  Index index() const {return idx_;}
  int first() const {return first_;}
  int last() const {return last_;}
  Node *parent() const {return parent_;}
  Node *child(int i) const {return child_[i];}
  hpx_addr_t sema() const {return sema_;}
  hpx_addr_t complete() const {return complete_;}

  void set_index(Index idx) {idx_ = idx;}
  void set_first(int first) {first_ = first;}
  void set_last(int last) {last_ = last;}
  void set_parent(Node *parent) {parent_ = parent;}
  void set_child(int i, Node *child) {child_[i] = child;} 
  void set_sema(hpx_addr_t sema) {sema_ = sema;}
  void set_complete(hpx_addr_t complete) {complete_ = complete;}
  
  void partition(Point *P, int *swap, int *bin, int *map, int threshold, 
                 double corner_x, double corner_y, double corner_z, 
                 double size);
  int n_descendants() const; 
  void compress(int *branch, int *tree, int parent, int &curr) const; 
  void extract(const int *branch, const int *tree, int n_nodes); 
  void destroy(bool allocated_in_array); 

private: 
  Index idx_; 
  int first_; 
  int last_; 
  Node *parent_; 
  Node *child_[8]; 
  hpx_addr_t sema_;
  hpx_addr_t complete_; 
};

#endif

