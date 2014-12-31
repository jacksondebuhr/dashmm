// =============================================================================
//  DASHMM 
//
//  Draft of tree construction for DASHMM
//  dashmm-tree.c
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


#include "dashmm-tree.h"

#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "libsync/sync.h"


//DEBUG
#include <stdio.h>


//These are used for some index calculations
// _offsets are the sums of the powers of 8
// _onlevels are the powers of 8
// _onlevel1d are the powers of 2
#define MAX_REFINE_LEVEL 8
const int _offsets[9]   = {0, 1, 9,  73,  585,  4681,  37449,  299593, 2396745};
const int _onlevels[8]  = {1, 8, 64, 512, 4096, 32768, 262144, 2097152};
const int _onlevel1d[8] = {1, 2, 4,  8,   16,   32,    64,     128};



//******* actions and action parameters *****************************


//NOTE: The naming scheme for the actions, the functions implementing the
// actions, and the parameters to those actions is pretty straightforward,
// if a bit long-winded. The action name has _action appended to the name 
// of the function implementing that action. The parameter type for that
// action has _params_t appended to the name of the function implementing
// the action.


hpx_action_t dashmm_tree_node_refine_action;

typedef struct {
  int refine_limit;
} dashmm_tree_node_refine_params_t;


hpx_action_t dashmm_tree_points_refine_action;

typedef struct {
  dashmm_volume_t vol;
  hpx_addr_t points;  //this is not strictly needed...
  int n_points;
} dashmm_tree_points_refine_params_t;

typedef struct {
  hpx_addr_t offsets[8];
  int n_points[8];
} dashmm_tree_points_refine_return_t;


hpx_action_t dashmm_tree_init_child_action;

typedef struct {
  dashmm_volume_t vol;
  hpx_addr_t parent;
  hpx_addr_t points;
  int n_points;
  int refine_limit;
} dashmm_tree_init_child_params_t;


hpx_action_t dashmm_parallel_record_spawn_action;

typedef struct {
  hpx_addr_t base_record;       //these update with the blocks
  int n_start;                  //the index of the record at base_record
  int n_records;                //...
  int64_t record_size;
  uint32_t block_size;
  hpx_action_t terminal_action;
  int branching_factor;
  int chunk_size;
  hpx_addr_t donelco;
  int64_t arg_size;
  char args[];
} dashmm_parallel_record_spawn_params_t;

typedef struct {
  hpx_addr_t donelco;
  int index;
  int64_t arg_size;
  char args[];
} dashmm_parallel_record_spawn_terminal_params_t;


hpx_action_t dashmm_tree_topnode_count_action;

typedef struct {
  dashmm_tree_t tree;
} dashmm_tree_topnode_count_params_t;


hpx_action_t dashmm_tree_topnode_init_action;

typedef struct {
  dashmm_tree_t tree;
} dashmm_tree_topnode_init_params_t;


hpx_action_t dashmm_tree_topnode_increment_count_action;

typedef struct {
  hpx_addr_t donelco;
  hpx_addr_t point;
} dashmm_tree_topnode_increment_count_params_t;


hpx_action_t dashmm_tree_topnode_points_alloc_action;

typedef struct {
  hpx_addr_t done_counting;
} dashmm_tree_topnode_points_alloc_params_t;


hpx_action_t dashmm_tree_topnode_start_refine_action;

typedef struct {
  int refine_limit;
} dashmm_tree_topnode_start_refine_params_t;


hpx_action_t dashmm_tree_destroy_node_action;

typedef struct {
  int top_depth;
  int this_depth;
} dashmm_tree_destroy_node_params_t;



//******* user interface functions **********************************

hpx_addr_t dashmm_tree_create(hpx_addr_t points,    //the points in GAS
                              int n_points,         //the number of points
                              int n_per_block,      //number of records per 
                                                    // block (for now?)
                              dashmm_volume_t vol,  //the volume (for now)
                              int refinement,       //the refinement limit
                              int top_depth         //the depth of the top
                              ) {
  //create the tree object
  dashmm_tree_t tree;
  tree.vol = vol;
  tree.refinement = refinement;
  tree.top_depth = top_depth;
  
  assert(top_depth < MAX_REFINE_LEVEL);
  tree.n_topnodes = _offsets[tree.top_depth+1];
  
  tree.topnodes = hpx_gas_global_alloc(tree.n_topnodes, 
                                       sizeof(dashmm_tree_node_t));
  assert(tree.topnodes != HPX_NULL);
  
  
  //put the tree object into GAS
  hpx_addr_t retval = hpx_gas_alloc(sizeof(dashmm_tree_t));
  assert(retval != HPX_NULL);
  
  hpx_addr_t tree_put_done = hpx_lco_future_new(0);
  assert(tree_put_done != HPX_NULL);
  
  hpx_gas_memput(retval, &tree, sizeof(tree), HPX_NULL, tree_put_done);
  
  
  //create the tree topnodes
  // this is a parallel spawn to fill in some values for the topnodes
  // the lowest topnodes will allocate a future for their first child
  // to act as the signal to start adding yourself.
  dashmm_tree_topnode_init_params_t initparms;
  memcpy(&initparms.tree, &tree, sizeof(dashmm_tree_t));
  
  hpx_addr_t topnode_done = dashmm_parallel_spawn(
                        tree.topnodes, tree.n_topnodes,
                        sizeof(dashmm_tree_node_t), sizeof(dashmm_tree_node_t),
                        dashmm_tree_topnode_init_action,
                        TREE_CREATE_BRANCHING_FACTOR, TREE_CREATE_CHUNK_SIZE,
                        &initparms, sizeof(initparms));
  hpx_lco_wait(topnode_done);
  hpx_lco_delete(topnode_done, HPX_NULL);
  //THIS IS A GLOBAL BARRIER - IS THERE A WAY TO REMOVE IT?
  
 
  //start the process of sending counts to the topnodes
  // using the parallel spawn
  dashmm_tree_topnode_count_params_t input;
  memcpy(&input.tree, &tree, sizeof(dashmm_tree_t));
    
  hpx_addr_t count_done = dashmm_parallel_spawn(
                points, n_points, 
                sizeof(dashmm_point_t), sizeof(dashmm_point_t) * n_per_block,
                dashmm_tree_topnode_count_action, 
                TREE_CREATE_BRANCHING_FACTOR, TREE_CREATE_CHUNK_SIZE,
                &input, sizeof(input));
  
  
  //start the process of allocating space for the points of the topnodes
  // using the parallel spawn - these will wait on count_done before doing
  // actual work
  dashmm_tree_topnode_points_alloc_params_t allocparms;
  allocparms.done_counting = count_done;
  
  hpx_addr_t finesttopnodes = 
                hpx_addr_add(tree.topnodes, 
                             sizeof(dashmm_tree_node_t) * _offsets[top_depth],
                             sizeof(dashmm_tree_node_t));
  int finestcount = _onlevels[top_depth];
  
  hpx_addr_t alloc_done = dashmm_parallel_spawn(
                finesttopnodes, finestcount,
                sizeof(dashmm_tree_node_t), sizeof(dashmm_tree_node_t),
                dashmm_tree_topnode_points_alloc_action,
                TREE_CREATE_BRANCHING_FACTOR, TREE_CREATE_CHUNK_SIZE,
                &allocparms, sizeof(allocparms));
  
  
  //start the process of refining the topnodes
  // this will use parallel spawn - this will wait on the LCO stored in the
  // second child before doing anything. 
  dashmm_tree_topnode_start_refine_params_t refparms;
  refparms.refine_limit = refinement;
  
  hpx_addr_t refine_done = dashmm_parallel_spawn(
                finesttopnodes, finestcount,
                sizeof(dashmm_tree_node_t), sizeof(dashmm_tree_node_t),
                dashmm_tree_topnode_start_refine_action,
                TREE_CREATE_BRANCHING_FACTOR, TREE_CREATE_CHUNK_SIZE,
                &refparms, sizeof(refparms));
   
  //cleanup
  hpx_lco_wait(alloc_done);
  hpx_lco_delete(alloc_done, HPX_NULL);
  
  hpx_lco_wait(refine_done);
  hpx_lco_delete(refine_done, HPX_NULL);
  
  hpx_lco_delete(count_done, HPX_NULL);

  
  //make sure the tree has been put into GAS
  hpx_lco_wait(tree_put_done);
  hpx_lco_delete(tree_put_done, HPX_NULL);
  
  return retval;
}


void dashmm_tree_destroy(hpx_addr_t tree_gas) {
  //get the tree from memory
  dashmm_tree_t tree;
  hpx_addr_t getdone = hpx_lco_future_new(0);
  assert(getdone != HPX_NULL);
  hpx_gas_memget(&tree, tree_gas, sizeof(tree), getdone);
  hpx_lco_wait(getdone);
  hpx_lco_delete(getdone, HPX_NULL);

  //This is simply an invocation of the tree destroy action
  dashmm_tree_destroy_node_params_t input;
  input.top_depth = tree.top_depth;
  input.this_depth = 0;
  hpx_call_sync(tree.topnodes, dashmm_tree_destroy_node_action, 
                &input, sizeof(input),
                NULL, 0);

  //used to make sure the tree is totally out of memory
  hpx_addr_t deldone = hpx_lco_and_new(2);

  //now that is done, we can get rid of the topnodes
  hpx_gas_free(tree.topnodes, deldone);
  
  //and the tree itself
  hpx_gas_free(tree_gas, deldone);
  
  hpx_lco_wait(deldone);
  hpx_lco_delete(deldone, HPX_NULL);
  
  return;
}


//******* utility functions *****************************************

void dashmm_tree_register_actions(void) {
  HPX_REGISTER_ACTION(&dashmm_tree_node_refine_action, 
                      dashmm_tree_node_refine);
  HPX_REGISTER_ACTION(&dashmm_tree_points_refine_action, 
                      dashmm_tree_points_refine);
  HPX_REGISTER_ACTION(&dashmm_tree_init_child_action, 
                      dashmm_tree_init_child);
  HPX_REGISTER_ACTION(&dashmm_parallel_record_spawn_action, 
                      dashmm_parallel_record_spawn);
  HPX_REGISTER_ACTION(&dashmm_tree_topnode_count_action, 
                      dashmm_tree_topnode_count);
  HPX_REGISTER_ACTION(&dashmm_tree_topnode_init_action, 
                      dashmm_tree_topnode_init);
  HPX_REGISTER_ACTION(&dashmm_tree_topnode_increment_count_action, 
                      dashmm_tree_topnode_increment_count);
  HPX_REGISTER_ACTION(&dashmm_tree_topnode_points_alloc_action, 
                      dashmm_tree_topnode_points_alloc);
  HPX_REGISTER_ACTION(&dashmm_tree_topnode_start_refine_action, 
                      dashmm_tree_topnode_start_refine);
  HPX_REGISTER_ACTION(&dashmm_tree_destroy_node_action, 
                      dashmm_tree_destroy_node);
}


int dashmm_volume_which_child(dashmm_volume_t vol, double *pos) {
  int result = 0;
  
  double center[3];
  center[0] = 0.5 * (vol.a[0] + vol.b[0]);
  center[1] = 0.5 * (vol.a[1] + vol.b[1]);
  center[2] = 0.5 * (vol.a[2] + vol.b[2]);
  
  if(pos[0] >= center[0]) {
    result += 1;
  }
  
  if(pos[1] >= center[1]) {
    result += 2;
  }
  
  if(pos[2] >= center[2]) {
    result += 4;
  }
  
  return result;
}


dashmm_volume_t dashmm_volume_of_child(dashmm_volume_t vol, int which) {
  dashmm_volume_t result;
  
  double center[3];
  center[0] = 0.5 * (vol.a[0] + vol.b[0]);
  center[1] = 0.5 * (vol.a[1] + vol.b[1]);
  center[2] = 0.5 * (vol.a[2] + vol.b[2]);
  
  if(which & 1) {
    result.a[0] = center[0];
    result.b[0] = vol.b[0];
  } else {
    result.a[0] = vol.a[0];
    result.b[0] = center[0];
  }
  
  if(which & 2) {
    result.a[1] = center[1];
    result.b[1] = vol.b[1];
  } else {
    result.a[1] = vol.a[1];
    result.b[1] = center[1];
  }
  
  if(which & 4) {
    result.a[2] = center[2];
    result.b[2] = vol.b[2];
  } else {
    result.a[2] = vol.a[2];
    result.b[2] = center[2];
  }
  
  return result;
}


void dashmm_invoke_tree_init_child_sync(hpx_addr_t node,
                                   dashmm_volume_t vol, hpx_addr_t parent,
                                   hpx_addr_t points, int n_points, 
                                   int refine_limit) {
  dashmm_tree_init_child_params_t input;
  input.vol = vol;
  input.parent = parent;
  input.points = points;
  input.n_points = n_points;
  input.refine_limit = refine_limit;
  
  hpx_call_sync(node, dashmm_tree_init_child_action, 
                &input, sizeof(input),
                NULL, 0);
}


int dashmm_tree_topnode_level(int index) {
  int retval = 0;
  
  for(int i = 0; i < MAX_REFINE_LEVEL; ++i) {
    if(index >= _offsets[i] && index < _offsets[i+1]) {
      retval = i;
      break;
    }
  }
  
  return retval;
}


void dashmm_tree_topnode_index_in_level(int index, int level, int *array) {
  //get the portion of index that is offset from the first at level
  int idx = index - _offsets[level];
  int onlevel = _onlevel1d[level];
  
  //now divide that index into chunks
  //NOTE: If we change the layout, this will change
  array[0] = idx % onlevel;
  idx -= array[0];
  idx /= onlevel;
  array[1] = idx % onlevel;
  idx -= array[1];
  idx /= onlevel;
  array[2] = idx;
}


dashmm_volume_t dashmm_tree_topnode_volume(dashmm_volume_t vol, 
                                           int level, int *index) {
  //vol here is the base volume
  double denom = 1.0 / ((double)_onlevel1d[level]);
  double factor[3];
  factor[0] = (vol.b[0] - vol.a[0]) * denom;
  factor[1] = (vol.b[1] - vol.a[1]) * denom;
  factor[2] = (vol.b[2] - vol.a[2]) * denom;
  
  dashmm_volume_t retval;
  retval.a[0] = vol.a[0] + factor[0] * index[0];
  retval.a[1] = vol.a[1] + factor[1] * index[1];
  retval.a[2] = vol.a[2] + factor[2] * index[2];
  retval.b[0] = retval.a[0] + factor[0];
  retval.b[1] = retval.a[1] + factor[1];
  retval.b[2] = retval.a[2] + factor[2];
  
  return retval;
}


hpx_addr_t dashmm_tree_topnode_address(hpx_addr_t base, int level, int *index) {
  //get the starting offset
  int offset = _offsets[level];
  
  //get the on level offset
  int onlevel = index[0] 
                + index[1] * _onlevel1d[level] 
                + index[2] * _onlevel1d[level] * _onlevel1d[level];
  
  offset += onlevel;
  
  //get the address
  hpx_addr_t retval = hpx_addr_add(base, 
                                   offset * sizeof(dashmm_tree_node_t),
                                   sizeof(dashmm_tree_node_t));
  
  return retval;
}


void dashmm_tree_topnode_from_point(dashmm_volume_t vol, int level,
                                    dashmm_point_t *point, int *index) {
  double denom = 1.0 / ((double)_onlevel1d[level]);
  double factor[3];
  factor[0] = (vol.b[0] - vol.a[0]) * denom;
  factor[1] = (vol.b[1] - vol.a[1]) * denom;
  factor[2] = (vol.b[2] - vol.a[2]) * denom;
  
  double num[3];
  num[0] = (point->pos[0] - vol.a[0]) / factor[0];
  num[1] = (point->pos[1] - vol.a[1]) / factor[1];
  num[2] = (point->pos[2] - vol.a[2]) / factor[2];
  
  index[0] = floor(num[0]);
  index[1] = floor(num[1]);
  index[2] = floor(num[2]);
  
  //DEBUG Assertions
  assert(index[0] < _onlevel1d[level]);
  assert(index[1] < _onlevel1d[level]);
  assert(index[2] < _onlevel1d[level]);
}


hpx_addr_t dashmm_parallel_spawn(
                        hpx_addr_t base_record,       //first record
                        int n_records,                //number of records
                        int64_t record_size,          //types match HPX-5
                        uint32_t block_size,          // hpx_addr_add
                        hpx_action_t terminal_action, //action to take
                        int branching_factor,         //tree branches 
                        int chunk_size,               //will not refine this
                        void *args,
                        int64_t arg_size) {
  //For the time being, return an and LCO as the termination detection
  hpx_addr_t retval = hpx_lco_and_new(n_records);
  assert(retval != HPX_NULL);
  
  //start the thing off
  size_t sz = sizeof(dashmm_parallel_record_spawn_params_t) + arg_size;
  dashmm_parallel_record_spawn_params_t *input = malloc(sz);
  input->base_record = base_record;
  input->n_start = 0;
  input->n_records = n_records;
  input->record_size = record_size;
  input->block_size = block_size;
  input->terminal_action = terminal_action;
  input->branching_factor = branching_factor;
  input->chunk_size = chunk_size;
  input->donelco = retval;
  input->arg_size = arg_size;
  if(args != NULL) {
    memcpy(input->args, args, arg_size);
  }
  
  hpx_call(base_record, dashmm_parallel_record_spawn_action,
           input, sz, HPX_NULL);
  
  //now that we are done with the input, we can rid ourselves of it
  free(input);
  
  //done
  return retval;
}


//******* actions ***************************************************

int dashmm_tree_node_refine(void *args) {
  hpx_addr_t node_gas = hpx_thread_current_target();
  dashmm_tree_node_t *node;
  if( !hpx_gas_try_pin(node_gas, (void **)&node) ) {
    return HPX_RESEND;
  }
  
  //interpret arguments
  dashmm_tree_node_refine_params_t *parms = args;
  
  //decide if we have too many particles
  if(node->n_points > parms->refine_limit) {
    //sort point sequence over the children
    dashmm_tree_points_refine_params_t input;
    input.vol = node->vol;
    input.points = node->points;
    input.n_points = node->n_points;
    
    dashmm_tree_points_refine_return_t result;
    
    hpx_call_sync(node->points, dashmm_tree_points_refine_action,
                  &input, sizeof(input),
                  &result, sizeof(result));
    
    //allocate memory for the children of the node, and count children
    int child_count = 0;
    int test_count = 0;
    for(int i = 0; i < 8; ++i) {
      if(result.n_points[i] > 0) {
        ++child_count;
        node->child[i] = hpx_gas_alloc(sizeof(dashmm_tree_node_t));
        assert(node->child[i] != HPX_NULL);
      }
      test_count += result.n_points[i];
    }
    assert(test_count == node->n_points);
    
    //create the children with their segments as arguments
    //spawn refinement at the child nodes
    //TODO assess if this synchronization is the best one, or if we want to 
    // relax it. For the moment this action will not return until all the 
    // children of this node have finished.
    
    hpx_addr_t alldone = hpx_lco_and_new(child_count);
    assert(alldone != HPX_NULL);
    
    for(int i = 0; i < 8; ++i) {
      if(result.n_points[i] > 0) {
        //spawn the action to create the child node
        dashmm_tree_init_child_params_t input;
        input.vol = dashmm_volume_of_child(node->vol, i);
        input.parent = node_gas;
        input.points = result.offsets[i];
        input.n_points = result.n_points[i];
        input.refine_limit = parms->refine_limit;
        
        hpx_call(node->child[i], dashmm_tree_init_child_action,
                 &input, sizeof(input), alldone);
      }
    }
    
    hpx_lco_wait(alldone);
    hpx_lco_delete(alldone, HPX_NULL);
    
  } else {
    //we are happy how we are
    //FOR NOW: nothing here, but in the future here is where we can start
    // the moment computation
    //FUTURE: We might want to create a local copy of the point data, so that
    // some operations can proceed quickly. Or not.
    
    //THIS branch is also executed when the node has no points. This is
    // relevant for topnodes that did not end up with any points.
  }
  
  //done
  hpx_gas_unpin(node_gas);
  return HPX_SUCCESS;
}


int dashmm_tree_points_refine(void *args) {
  hpx_addr_t points_gas = hpx_thread_current_target();
  dashmm_point_t *points;
  if( !hpx_gas_try_pin(points_gas, (void **)&points) ) {
    return HPX_RESEND;
  }
  
  dashmm_tree_points_refine_params_t *parms = args;
  
  //setup defaults on the result
  dashmm_tree_points_refine_return_t result;
  for(int i = 0; i < 8; ++i) {
    result.offsets[i] = points_gas;
    result.n_points[i] = 0;
  }
  
  //find the bin for each
  for(int i = 0; i < parms->n_points; ++i) {
    int sort = dashmm_volume_which_child(parms->vol, points[i].pos);
    points[i].sort = sort;
    result.n_points[sort] += 1;
  }
  
  //setup offsets for the different blocks
  int place[9] = {0};
  int total_offset = 0;
  for(int i = 0; i < 8; ++i) {
    place[i] = total_offset;
    result.offsets[i] = hpx_addr_add(points_gas, 
                                    total_offset * sizeof(dashmm_point_t),
                                    parms->n_points * sizeof(dashmm_point_t));
    total_offset += result.n_points[i];
  }
  place[9] = total_offset;
  
  //first give each particle an ordering
  for(int i = 0; i < parms->n_points; ++i) {
    int sort = points[i].sort;
    points[i].sort = place[sort];
    place[sort] += 1;
  }
  
  //TODO: we can probably do better here. This will continue to sort even
  // if the particles are in the same block. Perhaps there is a better way
  // to do this.
  // Though, as is, this sorting is stable. So if we want that property, then
  // we should keep this as it is now.
  
  //reorder the particles, this moves the particles at most once each
  dashmm_point_t source, save;
  int isource, isave, idest;
  
  for(int i = 0; i < parms->n_points; ++i) {
    if(points[i].sort != i) {
      memcpy(&source, &points[i], sizeof(source));
      isource = points[i].sort;
      
      idest = points[i].sort;
      
      do {
        memcpy(&save, &points[idest], sizeof(save));
        isave = points[idest].sort;
        
        memcpy(&points[idest], &source, sizeof(source));
        
        if(idest == i) {
          break;
        }
        
        memcpy(&source, &save, sizeof(source));
        isource = isave;
        
        idest = isource;
      } while(1);
    }
  }
  
  //all done
  hpx_gas_unpin(points_gas);
  
  //continue the result
  hpx_thread_continue(sizeof(result), &result);
}


//This can be used in two modes. If the refine_limit parameter is set in the
//input arguments, this will call the node refinement action synchronously.
//Otherwise, this will just set the values before finishing.
int dashmm_tree_init_child(void *args) {
  hpx_addr_t node_gas = hpx_thread_current_target();
  dashmm_tree_node_t *node;
  if( !hpx_gas_try_pin(node_gas, (void **)&node) ) {
    return HPX_RESEND;
  }
  
  dashmm_tree_init_child_params_t *parms = args;
  
  //set some things
  node->vol = parms->vol;
  node->parent = parms->parent;
  node->points = parms->points;
  node->n_points = parms->n_points;
  
  for(int i = 0; i < 8; ++i) {
    node->child[i] = HPX_NULL;
  }
  
  hpx_gas_unpin(node_gas);
  
  
  //if the refine_limit is nonzero, refine the node
  if(parms->refine_limit) {
    hpx_call_sync(node_gas, dashmm_tree_node_refine_action,
                  &parms->refine_limit, sizeof(parms->refine_limit),
                  NULL, 0);
  }
  
  
  //done
  return HPX_SUCCESS;
}


int dashmm_parallel_record_spawn(void *args) {
  //interpret params
  dashmm_parallel_record_spawn_params_t *parms = args;
  
  if(parms->n_records <= parms->chunk_size) {
    //set up the arguments...
    dashmm_parallel_record_spawn_terminal_params_t *input;
    size_t sz = sizeof(dashmm_parallel_record_spawn_terminal_params_t)
                    + parms->arg_size;
    input = malloc(sz);
    assert(input != NULL);
    
    if(parms->args != NULL) {
      memcpy(input->args, parms->args, parms->arg_size);
    }
    input->arg_size = parms->arg_size;
    input->donelco = parms->donelco;
  
    //just do the straight up loop here
    for(int irec = 0; irec < parms->n_records; ++irec) {
      //compute the target addess
      input->index = parms->n_start + irec;
      hpx_addr_t target = hpx_addr_add(parms->base_record, 
                                      irec * parms->record_size,
                                      parms->block_size);
      
      //invoke the action
      hpx_call(target, parms->terminal_action, input, sz, HPX_NULL);
    }
  } else {
    //break the remaining into pieces
    int n_chunks = parms->branching_factor;
    int ck_size = parms->n_records / n_chunks;
    if(parms->n_records % n_chunks) {
      ++ck_size;
    }
    
    //if we do not have enough records, we do something a bit different
    // and compute the best number of chunks to use
    if(ck_size < parms->chunk_size) {
      n_chunks = parms->n_records / parms->chunk_size;
      if(parms->n_records % parms->chunk_size) {
        ++n_chunks;
      }
      ck_size = parms->chunk_size;
    }
    
    hpx_addr_t oldbase = parms->base_record;
    int oldnstart = parms->n_start;
    int oldnrecs = parms->n_records;
    
    for(int i = 0; i < n_chunks; ++i) {
      int offset = i * ck_size;
      if(offset + ck_size > parms->n_records) {
        //this will only happen on the last iteration, so it is fine to 
        // modify ck_size
        ck_size = oldnrecs - i * ck_size;
      }
      
      //setup the arguments
      parms->base_record = hpx_addr_add(oldbase, 
                                        offset * parms->record_size,
                                        parms->block_size);
      parms->n_start = oldnstart + offset;
      parms->n_records = ck_size;
      
      //invoke the action
      hpx_call(parms->base_record, dashmm_parallel_record_spawn_action,
               parms, sizeof(*parms) + parms->arg_size, HPX_NULL);
    }
  }
  
  return HPX_SUCCESS;
}


int dashmm_tree_topnode_count(void *args) {
  hpx_addr_t point_gas = hpx_thread_current_target();
  dashmm_point_t *point;
  if( !hpx_gas_try_pin(point_gas, (void **)&point) ) {
    return HPX_RESEND;
  }
  
  dashmm_parallel_record_spawn_terminal_params_t *wrap = args;
  dashmm_tree_topnode_count_params_t *parms = 
                (dashmm_tree_topnode_count_params_t *)wrap->args;
  
  //decide which topnode I am in
  // This takes the point data, the volume
  int pidx[3];
  dashmm_tree_topnode_from_point(parms->tree.vol, parms->tree.top_depth,
                                            point, pidx);
  hpx_addr_t target = dashmm_tree_topnode_address(parms->tree.topnodes,
                                    parms->tree.top_depth, pidx);
  
  //inform the topnode of the point
  dashmm_tree_topnode_increment_count_params_t input;
  input.donelco = wrap->donelco;
  input.point = point_gas;
  
  hpx_call(target, dashmm_tree_topnode_increment_count_action, 
                           &input, sizeof(input), HPX_NULL);
  
  hpx_gas_unpin(point_gas);
  
  return HPX_SUCCESS;
}


int dashmm_tree_topnode_increment_count(void *args) {
  hpx_addr_t node_gas = hpx_thread_current_target();
  dashmm_tree_node_t *node;
  if( !hpx_gas_try_pin(node_gas, (void **)&node) ) {
    return HPX_RESEND;
  }
  
  //increment atomically
  sync_fadd(&(node->n_points), 1, SYNC_ACQ_REL);
  
  //interpret parameters
  dashmm_tree_topnode_increment_count_params_t *parms = args;
  
  //set the counting lco
  hpx_lco_and_set(parms->donelco, HPX_NULL);
  
  //wait for the first lco to be set
  hpx_lco_wait(node->child[0]);
  
  //atomically get place in index array
  int dest_index = sync_fadd(&(node->n_arrived), 1, SYNC_ACQ_REL);
  
  //move memory about
  hpx_addr_t target = hpx_addr_add(node->points, 
                                   dest_index * sizeof(dashmm_point_t),
                                   node->n_points * sizeof(dashmm_point_t));
  hpx_addr_t movedone = hpx_lco_future_new(0);
  assert(movedone != HPX_NULL);
  
  hpx_gas_memcpy(target, parms->point, sizeof(dashmm_point_t), movedone);
  
  hpx_lco_wait(movedone);
  hpx_lco_delete(movedone, HPX_NULL);
  
  //set an input on the and gate
  hpx_lco_and_set(node->child[1], HPX_NULL);
  
  //done
  hpx_gas_unpin(node_gas);
  return HPX_SUCCESS;
}


int dashmm_tree_topnode_init(void *args) {
  hpx_addr_t node_gas = hpx_thread_current_target();
  dashmm_tree_node_t *node;
  if( !hpx_gas_try_pin(node_gas, (void **)&node) ) {
    return HPX_RESEND;
  }
  
  //interpret the parameters
  dashmm_parallel_record_spawn_terminal_params_t *wrap = args;
  dashmm_tree_topnode_init_params_t *parms = 
            (dashmm_tree_topnode_init_params_t *)wrap->args;
  
  //fill in the various needed parameters and so on
  int level = dashmm_tree_topnode_level(wrap->index);
  int index[3];
  dashmm_tree_topnode_index_in_level(wrap->index, level, index);
  node->vol = dashmm_tree_topnode_volume(parms->tree.vol, level, index);
  
  //UM: What is going on here. This sets the parent to itself...
  node->parent = dashmm_tree_topnode_address(parms->tree.topnodes, 
                                                        level, index);
  assert(node->parent != node_gas);
  node->n_points = 0;
  node->n_arrived = 0;
  
  node->points = HPX_NULL;

  if(level < parms->tree.top_depth) {
    int clevel = level + 1;
    int cindex[3];
    cindex[0] = 2 * index[0];
    cindex[1] = 2 * index[1];
    cindex[2] = 2 * index[2];
    
    node->child[0] = dashmm_tree_topnode_address(parms->tree.topnodes,
                                                    clevel, cindex);
    cindex[0] += 1;
    node->child[1] = dashmm_tree_topnode_address(parms->tree.topnodes,
                                                    clevel, cindex);
    cindex[0] -= 1;
    cindex[1] += 1;
    node->child[2] = dashmm_tree_topnode_address(parms->tree.topnodes,
                                                    clevel, cindex);
    cindex[0] += 1;
    node->child[3] = dashmm_tree_topnode_address(parms->tree.topnodes,
                                                    clevel, cindex);
    cindex[0] -= 1;
    cindex[1] -= 1;
    cindex[2] += 1;
    node->child[4] = dashmm_tree_topnode_address(parms->tree.topnodes,
                                                    clevel, cindex);
    cindex[0] += 1;
    node->child[5] = dashmm_tree_topnode_address(parms->tree.topnodes,
                                                    clevel, cindex);
    cindex[0] -= 1;
    cindex[1] += 1;
    node->child[6] = dashmm_tree_topnode_address(parms->tree.topnodes,
                                                    clevel, cindex);
    cindex[0] += 1;
    node->child[7] = dashmm_tree_topnode_address(parms->tree.topnodes,
                                                    clevel, cindex);
  } else {
    //the first is used to signal that the point memory has been allocated
    node->child[0] = hpx_lco_future_new(0);
    assert(node->child[0] != HPX_NULL);
    
    node->child[1] = HPX_NULL;
    node->child[2] = HPX_NULL;
    node->child[3] = HPX_NULL;
    node->child[4] = HPX_NULL;
    node->child[5] = HPX_NULL;
    node->child[6] = HPX_NULL;
    node->child[7] = HPX_NULL;
  }
  
  //set the done LCO
  hpx_lco_and_set(wrap->donelco, HPX_NULL);
  
  //unpin and done
  hpx_gas_unpin(node_gas);
  
  return HPX_SUCCESS;
}


int dashmm_tree_topnode_points_alloc(void *args) {
  hpx_addr_t node_gas = hpx_thread_current_target();
  dashmm_tree_node_t *node;
  if( !hpx_gas_try_pin(node_gas, (void **)&node) ) {
    return HPX_RESEND;
  } 
  
  //interpret parameters 
  dashmm_parallel_record_spawn_terminal_params_t *wrap = args;
  dashmm_tree_topnode_points_alloc_params_t *parms = 
            (dashmm_tree_topnode_points_alloc_params_t *)wrap->args;
  
  //wait for the counting to be done
  hpx_lco_wait(parms->done_counting);
  
  if(node->n_points) {
    //allocate the and gate
    node->child[1] = hpx_lco_and_new(node->n_points);
    assert(node->child[1] != HPX_NULL);
  
    //allocate space
    node->points = hpx_gas_alloc(sizeof(dashmm_point_t) * node->n_points);
    assert(node->points != HPX_NULL);
  
    //set LCO to indicate that it is ready
    hpx_lco_set(node->child[0], 0, NULL, HPX_NULL, HPX_NULL);
  }
  
  //set the completion LCO as well
  hpx_lco_and_set(wrap->donelco, HPX_NULL);
  
  hpx_gas_unpin(node_gas);
  
  return HPX_SUCCESS;
}


int dashmm_tree_topnode_start_refine(void *args) {
  hpx_addr_t node_gas = hpx_thread_current_target();
  dashmm_tree_node_t *node;
  if( !hpx_gas_try_pin(node_gas, (void **)&node) ) {
    return HPX_RESEND;
  }
  
  //interpret params
  dashmm_parallel_record_spawn_terminal_params_t *wrap = args;
  dashmm_tree_topnode_start_refine_params_t *parms = 
            (dashmm_tree_topnode_start_refine_params_t *)wrap->args;
  
  //wait for all to have arrived
  hpx_lco_wait(node->child[0]);
  
  int refine = 0;
  if(node->n_points) {
    refine = 1;
    hpx_lco_wait(node->child[1]);
    hpx_lco_delete(node->child[1], HPX_NULL);
  }
  
  hpx_lco_delete(node->child[0], HPX_NULL);
  
  
  //done with the node data, so unpin
  hpx_gas_unpin(node_gas);
  
  if(refine) {
    //begin the refinement process
    hpx_call_sync(node_gas, dashmm_tree_node_refine_action,
                &parms->refine_limit, sizeof(parms->refine_limit),
                NULL, 0);
  }
  
  //set the done lco
  hpx_lco_and_set(wrap->donelco, HPX_NULL);
  
  return HPX_SUCCESS;
}


//NOTE: Something goes very wrong here... Not sure what it is...
// So don't call this routine it seems.
int dashmm_tree_destroy_node(void *args) {
  hpx_addr_t node_gas = hpx_thread_current_target();
  dashmm_tree_node_t *node;
  if( !hpx_gas_try_pin(node_gas, (void **)&node) ) {
    return HPX_RESEND;
  }
  
  //interpret params
  dashmm_tree_destroy_node_params_t *parms = args;
  
  //spawn at children
  dashmm_tree_destroy_node_params_t input;
  input.top_depth = parms->top_depth;
  input.this_depth = parms->this_depth + 1;
  
  //count children
  int ccount = 0;
  for(int i = 0; i < 8; ++i) {
    if(node->child[i] != HPX_NULL) {
      ++ccount;
    }
  }
  
  hpx_addr_t cdone = hpx_lco_and_new(ccount);
  assert(cdone != HPX_NULL);
  
  for(int i = 0; i < 8; ++i) {
    if(node->child[i] != HPX_NULL) {
      hpx_call(node->child[i], dashmm_tree_destroy_node_action, 
               &input, sizeof(input), cdone);
    }
  }
  
  hpx_lco_wait(cdone);
  hpx_lco_delete(cdone, HPX_NULL);
  
  //free the children
  if(parms->this_depth > parms->top_depth) {
    for(int i = 0; i < 8; ++i) {
      if(node->child[i] != HPX_NULL) {
        hpx_gas_free(node->child[i], HPX_NULL);
      }
    }
  }
  
  //if this is the bottom of the topnodes, free the points
  if(parms->this_depth == parms->top_depth) {
    hpx_gas_free(node->points, HPX_NULL);
  }
  
  //done
  hpx_gas_unpin(node_gas);
  
  return HPX_SUCCESS;
}


