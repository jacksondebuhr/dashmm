#ifndef __TREE_H__
#define __TREE_H__

#include <vector>
#include "hpx/hpx.h"

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
  Node(Index idx = Index{}): idx_{idx}, first_{0}, last_{0}, parent_{nullptr} {
    for (int i = 0; i < 8; ++i) 
      child_[i] = nullptr; 
    sema_ = hpx_lco_sema_new(1); 
  }

  Node(Index idx, int first, int last, Node *parent) : 
    idx_{idx}, first_{first}, last_{last}, parent_{parent}, is_leaf_{false}
  {
    for (int i = 0; i < 8; ++i) 
      child_[i] = nullptr; 
    sema_ = hpx_lco_sema_new(1); 
  }

  ~Node() {
    for (int i = 0; i < 8; ++i) {
      if (child_[i]) 
        delete child_[i];
    }

    hpx_lco_delete_sync(sema_);
  }

  Index index() const {return idx_;}
  int first() const {return first_;}
  int last() const {return last_;}
  Node *parent() const {return parent_;}
  Node *child(int i) const {return child_[i];}
  
private: 
  Index idx_; 
  int first_; 
  int last_; 
  Node *parent_; 
  Node *child_[8]; 
  std::vector<const Node *> list1_; 
  std::vector<const Node *> list2_; 
  std::vector<const Node *> list3_; 
  std::vector<const Node *> list4_; 
  hpx_addr_t sema_;
} Node; 

extern hpx_action_t allocate_points_action; 
extern hpx_action_t set_points_action; 
extern hpx_action_t delete_points_action; 
extern hpx_action_t partition_points_action;   
extern hpx_action_t domain_geometry_init_action; 
extern hpx_action_t domain_geometry_op_action; 
extern hpx_action_t domain_geometry_predicate_action; 
extern hpx_action_t set_domain_geometry_action; 
extern hpx_action_t coarse_grid_count_init_action; 
extern hpx_action_t coarse_grid_count_op_action; 
extern hpx_action_t coarse_grid_count_predicate_action; 
extern hpx_action_t init_coarse_partition_action; 
extern hpx_action_t coarse_partition_action; 
extern hpx_action_t exchange_coarse_grid_count_action; 


int *distribute_points(int num_ranks, const int *global, int len); 



extern hpx_action_t delete_tree_action; 
extern hpx_action_t test_action; 
#endif 
