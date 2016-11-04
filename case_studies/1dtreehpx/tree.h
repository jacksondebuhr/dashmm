// ============================================================================
//  High Performance ParalleX Library (hpx-apps)
//
//  Copyright (c) 2013-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license.  See the COPYING file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// ============================================================================

#ifndef __1DTREE_TUTORIAL_TREE_H__
#define __1DTREE_TUTORIAL_TREE_H__


#include "hpx/hpx.h"


//void partition_node_sync(hpx_addr_t node, hpx_addr_t parts, int n_parts,
//                         int n_partition);

//void spawn_potential_computation(hpx_addr_t root, hpx_addr_t sync,
//                                 double theta);


//struct Moment {
//  double mtot;
//  double xcom;
//  double Q00;
//
//  Moment() : mtot{0.0}, xcom{0.0}, Q00{0.0} { }
//};


struct Particle {
  double pos;
  double mass;
  double phi;
};



struct Node {
 public:
  Node(double l, double h)
    : left{nullptr}, right{nullptr},
      low(l), high(h), parts{nullptr}, count{0} { }
  ~Node();


  void partition(Particle *parts, int n_parts, int thresh,
                 hpx_addr_t sync = HPX_NULL);


  // Action handlers
  static void register_actions();
  static int partition_handler(Node *node, Particle *parts, int n_parts,
                               int thresh);

  // Action ids
  static hpx_action_t partition_;

  // The data
  Node *left;
  Node *right;
  double low;
  double high;
  Particle *parts;
  int count;
};


#endif // __1DTREE_TUTORIAL_TREE_H__
