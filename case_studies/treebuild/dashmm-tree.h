// =============================================================================
//  DASHMM 
//
//  Draft of tree construction for DASHMM
//  dashmm-tree.h
//
//  Copyright (c) 2014, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license.  See the LICENSE file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
//
//  Authors:
//    Jackson DeBuhr, Indiana University <jdebuhr [at] indiana.edu>
// =============================================================================

#ifndef __DASHMM_TREE_H__
#define __DASHMM_TREE_H__


#include "hpx/hpx.h"


//TODO pick these in some kind of intelligent fashion
// These would need to be tuned based on actual performance I think.
#define TREE_CREATE_BRANCHING_FACTOR 2
#define TREE_CREATE_CHUNK_SIZE 10


//In the future, this is more likely going to be part of a dashmm_types.h
// or similar. Though, if it is not used elsewhere, then it could stay here.
typedef struct {
  double a[3];
  double b[3];
} dashmm_volume_t;


typedef struct {
  double pos[3];
  //eventually the charges or the computed forces - 
  //  probably a union so we can select based on context
  int sort;
} dashmm_point_t;


typedef struct {
  dashmm_volume_t vol;  //the volume of the node
  
  hpx_addr_t parent;    //parent of this node
  hpx_addr_t child[8];  //children of this node
  
  hpx_addr_t points;    //address of the first point represented by this node
  uint32_t n_points;    //number of points represented by this node
  uint32_t n_arrived;   //number of points that have arrived
  
  //eventually the moments and so on...
} dashmm_tree_node_t;


typedef struct {
  dashmm_volume_t vol;  //the volume of the root
  
  int refinement;       //the node refinement limit
  
  int top_depth;        //the depth of the top of the tree
  hpx_addr_t topnodes;  //the topnodes of the tree
  int n_topnodes;       //the number of topnodes
} dashmm_tree_t;


//NOTE: for the time being, we prototype everything here. Eventually, we would
// only expose those routines needed from outside.

//user interface
hpx_addr_t dashmm_tree_create(hpx_addr_t points,
                              int n_points,
                              int n_per_block,
                              dashmm_volume_t vol,
                              int refinement,
                              int top_depth);


//utility functions
void dashmm_tree_register_actions(void);
int dashmm_volume_which_child(dashmm_volume_t vol, double *pos);
dashmm_volume_t dashmm_volume_of_child(dashmm_volume_t vol, int which);
void dashmm_invoke_tree_init_child_sync(hpx_addr_t node,
                                   dashmm_volume_t vol, hpx_addr_t parent,
                                   hpx_addr_t points, int n_points, 
                                   int refine_limit);
int dashmm_tree_topnode_level(int index);
void dashmm_tree_topnode_index_in_level(int index, int level, int *array);
dashmm_volume_t dashmm_tree_topnode_volume(dashmm_volume_t vol, 
                                           int level, int *index);
hpx_addr_t dashmm_tree_topnode_address(hpx_addr_t base, int level, int *index);
void dashmm_tree_topnode_from_point(dashmm_volume_t vol, int level,
                                    dashmm_point_t *point, int *index);


//NOTE: This assumes that the terminal actions all take the same arguments
hpx_addr_t dashmm_parallel_spawn(
                        hpx_addr_t base_record,       //first record
                        int n_records,                //number of records
                        int64_t record_size,          //types match HPX-5
                        uint32_t block_size,          // hpx_addr_add
                        hpx_action_t terminal_action, //action to take
                        int branching_factor,         //tree branches 
                        int chunk_size,               //will not refine this
                        void *args,
                        int64_t arg_size);
                                        


//actions
int dashmm_tree_node_refine(void *args);
int dashmm_tree_points_refine(void *args);
int dashmm_tree_init_child(void *args);
int dashmm_parallel_record_spawn(void *args);
int dashmm_tree_topnode_count(void *args);
int dashmm_tree_topnode_init(void *args);
int dashmm_tree_topnode_increment_count(void *args);
int dashmm_tree_topnode_points_alloc(void *args);
int dashmm_tree_topnode_start_refine(void *args);

#endif
