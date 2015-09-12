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


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "libsync/sync.h"

//
// Module variables
//


//These are used for some index calculations
// _offsets are the sums of the powers of 8 - these give the index of the first
//     topnode at a given refinement level; or they give the total number of
//     topnodes if you use the next value. Ex: the first level 2 topnode is at
//     index 9, and there would be 73 total topnodes.
//
// _onlevels are the powers of 8 - these give the number of topnodes at the 
//     level. The cumulative sum of these values is in _offsets
//
// _onlevel1d are the powers of 2 - these are the number of topnodes at each 
//     level in any given direction. Ex: there are 4 divisions for two levels
//     of refinement
#define MAX_REFINE_LEVEL 8
const int _offsets[9]   = {0, 1, 9,  73,  585,  4681,  37449,  299593, 2396745};
const int _onlevels[8]  = {1, 8, 64, 512, 4096, 32768, 262144, 2097152};
const int _onlevel1d[8] = {1, 2, 4,  8,   16,   32,    64,     128};



//
// Action and action parameters
//


hpx_action_t dashmm_block_spawn_action;

typedef struct {
  hpx_addr_t base;
  int block_size;
  int block_count;
  int first_block;
  int last_block;
  int result_size;
  void (*reduce)(void *a, void *b);
  hpx_action_t terminal_action;
  int payload_size;
  char payload[];
} dashmm_block_spawn_params_t;

typedef struct {
  int block_size;
  int block_number;
  int payload_size;
  char payload[];
} dashmm_block_spawn_terminal_params_t;

hpx_action_t dashmm_point_volume_action;


hpx_action_t dashmm_tree_node_init_action;

typedef struct {
  dashmm_tree_t tree;
  int point_block_count;
} dashmm_tree_node_init_params_t;


hpx_action_t dashmm_tree_points_sort_action;


hpx_action_t dashmm_tree_inform_topnode_action;


hpx_action_t dashmm_tree_node_alloc_action;

typedef struct {
  hpx_addr_t sortdone;
  int refine_limit;
} dashmm_tree_node_alloc_params_t;


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
#ifdef DEBUGID
  uint64_t ID;
#endif
} dashmm_tree_init_child_params_t;


//
// User Interface
//


hpx_addr_t dashmm_tree_create(hpx_addr_t points,
                              int n_points,
                              int n_per_block,
                              int refinement,
                              int top_depth) {
  hpx_time_t vol_start = hpx_time_now();

  //first get the volume
  int n_blocks = n_points / n_per_block;
  //For now, we assume that the points fit into the blocks evenly...
  assert( (n_points % n_per_block) == 0 );
  hpx_addr_t volume = dashmm_parallel_block_spawn(points,
                            n_per_block * sizeof(dashmm_point_t),
                            n_blocks,
                            dashmm_point_volume_action,
                            0,
                            NULL,
                            sizeof(dashmm_volume_t),
                            dashmm_volume_reducer);

  dashmm_volume_t domain;
  hpx_lco_get(volume, sizeof(domain), &domain);
  
  //cubify the volume - leaving a bit of space around the outermost points
  dashmm_volume_cubify(&domain, 1.005);
  
  hpx_time_t vol_end = hpx_time_now();
  
  //Create the tree data
  dashmm_tree_t tree;
  tree.vol = domain;
  tree.refinement = refinement;
  tree.top_depth = top_depth;
  
  assert(top_depth < MAX_REFINE_LEVEL);
  tree.n_topnodes = _offsets[tree.top_depth + 1];
  
  tree.topnodes = hpx_gas_global_alloc(tree.n_topnodes, 
                                       sizeof(dashmm_tree_node_t));
  assert(tree.topnodes != HPX_NULL);
  
  //put it into GAS
  hpx_addr_t retval = hpx_gas_alloc(sizeof(dashmm_tree_t));
  assert(retval != HPX_NULL);
  
  hpx_addr_t tree_put_done = hpx_lco_future_new(0);
  assert(tree_put_done != HPX_NULL);
  
  hpx_gas_memput(retval, &tree, sizeof(tree), HPX_NULL, tree_put_done);
  
  
  //block spawn to initialize the topnodes
  hpx_time_t init_start = hpx_time_now();
  
  dashmm_tree_node_init_params_t initinput;
  initinput.tree = tree;
  initinput.point_block_count = n_blocks;
  hpx_addr_t initdone = dashmm_parallel_block_spawn(tree.topnodes,
                            sizeof(dashmm_tree_node_t),
                            tree.n_topnodes,
                            dashmm_tree_node_init_action,
                            sizeof(initinput),
                            &initinput,
                            0,
                            NULL);
  hpx_lco_wait(initdone);
  hpx_lco_delete(initdone, HPX_NULL);
  
  hpx_time_t init_end = hpx_time_now();
  
  //block spawn to begin sorting the point data
  hpx_time_t sort_start = hpx_time_now();
  hpx_addr_t sortdone = dashmm_parallel_block_spawn(points, 
                            n_per_block * sizeof(dashmm_point_t),
                            n_blocks,
                            dashmm_tree_points_sort_action,
                            sizeof(tree),
                            &tree,
                            0,
                            NULL);
  
  //block spawn to begin the allocation actions for the finest topnodes
  hpx_time_t refine_start = hpx_time_now();
  hpx_addr_t first_finest = hpx_addr_add(tree.topnodes, 
                    sizeof(dashmm_tree_node_t) * _offsets[tree.top_depth],
                    sizeof(dashmm_tree_node_t));
  int n_finest = _onlevels[tree.top_depth];
  dashmm_tree_node_alloc_params_t allocinput;
  allocinput.sortdone = sortdone;
  allocinput.refine_limit = tree.refinement;
  hpx_addr_t refinedone = dashmm_parallel_block_spawn(first_finest,
                            sizeof(dashmm_tree_node_t),
                            n_finest,
                            dashmm_tree_node_alloc_action,
                            sizeof(allocinput),
                            &allocinput,
                            0,
                            NULL);
  
  //wait on those two things - really only the second is needed, as it cannot
  // finish without the first.
  hpx_lco_wait(sortdone);
  hpx_time_t sort_end = hpx_time_now();
  
  hpx_lco_wait(refinedone);
  hpx_time_t refine_end = hpx_time_now();
  
  hpx_lco_delete(sortdone, HPX_NULL);
  hpx_lco_delete(refinedone, HPX_NULL);
  
  //some final work
  hpx_lco_wait(tree_put_done);
  hpx_lco_delete(tree_put_done, HPX_NULL);
  
  double dt_volume = hpx_time_diff_ms(vol_start, vol_end);
  double dt_init = hpx_time_diff_ms(init_start, init_end);
  double dt_sort = hpx_time_diff_ms(sort_start, sort_end);
  double dt_refine = hpx_time_diff_ms(refine_start, refine_end);
  double dt_deltasortrefine = hpx_time_diff_ms(sort_end, refine_end);
  
  fprintf(stdout, "%lg %lg %lg %lg %lg ", 
	   dt_volume, dt_init, dt_sort, dt_refine, dt_deltasortrefine);
  
  return retval;
}
                              

void dashmm_tree_destroy(hpx_addr_t tree_gas) {
  //TODO
}



//
// Utility functions
//


void dashmm_tree_register_actions(void) {
  HPX_REGISTER_ACTION(dashmm_block_spawn, &dashmm_block_spawn_action);
  HPX_REGISTER_ACTION(dashmm_point_volume, &dashmm_point_volume_action);
  HPX_REGISTER_ACTION(dashmm_tree_node_init, &dashmm_tree_node_init_action);
  HPX_REGISTER_ACTION(dashmm_tree_points_sort, &dashmm_tree_points_sort_action);
  HPX_REGISTER_ACTION(dashmm_tree_inform_topnode, 
                         &dashmm_tree_inform_topnode_action);
  HPX_REGISTER_ACTION(dashmm_tree_node_alloc, &dashmm_tree_node_alloc_action);
  HPX_REGISTER_ACTION(dashmm_tree_node_refine, &dashmm_tree_node_refine_action);
  HPX_REGISTER_ACTION(dashmm_tree_points_refine, 
                         &dashmm_tree_points_refine_action);
  HPX_REGISTER_ACTION(dashmm_tree_init_child, &dashmm_tree_init_child_action);
}


int dashmm_volume_which_child(dashmm_volume_t vol, double *pos) {
  int result = 0;
  
  double center[3];
  center[0] = 0.5 * (vol.a[0] + vol.b[0]);
  center[1] = 0.5 * (vol.a[1] + vol.b[1]);
  center[2] = 0.5 * (vol.a[2] + vol.b[2]);
  
  if (pos[0] >= center[0]) {
    result += 1;
  }
  
  if (pos[1] >= center[1]) {
    result += 2;
  }
  
  if (pos[2] >= center[2]) {
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
  
  if (which & 1) {
    result.a[0] = center[0];
    result.b[0] = vol.b[0];
  } else {
    result.a[0] = vol.a[0];
    result.b[0] = center[0];
  }
  
  if (which & 2) {
    result.a[1] = center[1];
    result.b[1] = vol.b[1];
  } else {
    result.a[1] = vol.a[1];
    result.b[1] = center[1];
  }
  
  if (which & 4) {
    result.a[2] = center[2];
    result.b[2] = vol.b[2];
  } else {
    result.a[2] = vol.a[2];
    result.b[2] = center[2];
  }
  
  return result;
}


void dashmm_volume_reducer(void *a, void *b) {
  dashmm_volume_t *save = a;
  dashmm_volume_t *mod = b;
  
  for (int dir = 0; dir < 3; ++dir) {
    if (save->a[dir] > mod->a[dir]) {
      save->a[dir] = mod->a[dir];
    } //no else, as save already has the right value
    
    if (save->b[dir] < mod->b[dir]) {
      save->b[dir] = mod->b[dir];
    } //no else, as save already has the right value
  }
}


void dashmm_volume_cubify(dashmm_volume_t *vol, double expand) {
  double diff[3];
  diff[0] = vol->b[0] - vol->a[0];
  diff[1] = vol->b[1] - vol->a[1];
  diff[2] = vol->b[2] - vol->a[2];
  
  double size = diff[0];
  if (diff[1] > size) {
    size = diff[1];
  }
  if (diff[2] > size) {
    size = diff[2];
  }
  size *= 0.5 * expand;
  
  double center[3];
  center[0] = 0.5 * (vol->a[0] + vol->b[0]);
  center[1] = 0.5 * (vol->a[1] + vol->b[1]);
  center[2] = 0.5 * (vol->a[2] + vol->b[2]);
  
  vol->a[0] = center[0] - size;
  vol->a[1] = center[1] - size;
  vol->a[2] = center[2] - size;
  vol->b[0] = center[0] + size;
  vol->b[1] = center[1] + size;
  vol->b[2] = center[2] + size;
}


int dashmm_tree_topnode_level(int index) {
  int retval = 0;
  
  for (int i = 0; i < MAX_REFINE_LEVEL; ++i) {
    if (index >= _offsets[i] && index < _offsets[i+1]) {
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


int dashmm_tree_topnode_offset_on_level(int level, int *array) {
  int retval = array[0];
  retval += array[1] * _onlevel1d[level];
  retval += array[2] * _onlevel1d[level] * _onlevel1d[level];
  return retval;
}


int dashmm_tree_topnode_offset_from_index(int level, int *array) {
  int retval = _offsets[level] 
                + dashmm_tree_topnode_offset_on_level(level, array);
  return retval;
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
  int offset = dashmm_tree_topnode_offset_from_index(level, index);
  
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


#ifdef DEBUGID
uint64_t dashmm_tree_node_id(int *index, int level) {
  int on_level_1d = 1 << level;
  
  int idonlevel = index[0] 
                    + on_level_1d * index[1] 
                    + on_level_1d * on_level_1d * index[2];
  
  uint64_t retval = 1000 * idonlevel + level;
  
  return retval;
}


uint64_t dashmm_tree_node_index_for_child(uint64_t id, int which) {
  uint64_t level = id % 1000;
  int on_level_1d = 1 << level;

  uint64_t work = id - level;
  work /= 1000;
  uint64_t idx = work % on_level_1d;
  work -= idx;
  work /= on_level_1d;
  uint64_t idy = work % on_level_1d;
  work -= idy;
  work /= on_level_1d;
  uint64_t idz = work;
  
  //move to next level
  idx *= 2;
  idy *= 2;
  idz *= 2;
  
  if (which & 1) {
    idx += 1;
  }
  if (which & 2) {
    idy += 1;
  }
  if (which & 4) {
    idz += 1;
  }
  
  uint64_t retval = idx 
                    + idy * on_level_1d 
                    + idz * on_level_1d * on_level_1d;
  retval *= 1000;
  retval += level + 1;
  
  return retval;
}
#endif


//NOTE: the incoming particles must all be labeled with a sort value that 
// gives the bin to which they are assigned. This means the entries in bins
// must match that number of points with the given bin marked, or this will
// behave badly.
void dashmm_point_bin_sort(dashmm_point_t *points, int n_points,
                        int *bins, int n_bins) {
  //int *bin_offsets = calloc(sizeof(int), n_bins);
  int *bin_offsets = malloc(sizeof(int) * n_bins);
  assert(bin_offsets != NULL);
  for (int ibin = 0; ibin < n_bins; ++ibin) {
    bin_offsets[ibin] = 0;
  }
  for (int ibin = 1; ibin < n_bins; ++ibin) {
    bin_offsets[ibin] = bin_offsets[ibin - 1] + bins[ibin - 1];
  }
  //consistency check
  assert( (bin_offsets[n_bins - 1] + bins[n_bins - 1]) == n_points );
  
  for (int ipt = 0; ipt < n_points; ++ipt) {
    int mybin = points[ipt].sort;
    points[ipt].sort = bin_offsets[mybin];
    bin_offsets[mybin] += 1;
  }
  
  free(bin_offsets);

  dashmm_point_t source, save;
  int isource, isave, idest;
  
  for(int i = 0; i < n_points; ++i) {
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
}


void dashmm_continue_cleanup(void *env) {
  free(env);
}


hpx_addr_t dashmm_parallel_block_spawn(hpx_addr_t base,
                                       int block_size,
                                       int block_count,
                                       hpx_action_t terminal_action,
                                       int payload_size,
                                       void *payload,
                                       int result_size,
                                       void (*reduce)(void *a, void *b)
                                       ) {
  hpx_addr_t retval = hpx_lco_future_new(result_size);
  assert(retval != HPX_NULL);
    
  if (block_count > 1) {
    size_t input_size = sizeof(dashmm_block_spawn_params_t) + payload_size;
    dashmm_block_spawn_params_t *input = malloc(input_size);
    assert(input != NULL);
    input->base = base;
    input->block_size = block_size;
    input->block_count = block_count;
    input->first_block = 0;
    input->last_block = block_count - 1;
    input->result_size = result_size;
    input->reduce = reduce;
    input->terminal_action = terminal_action;
    input->payload_size = payload_size;
    memcpy(input->payload, payload, payload_size);
  
    hpx_call(base, dashmm_block_spawn_action, input, input_size, retval);
    
    free(input);
  } else {
    //we handle the special case of only one block here
    size_t input_size = sizeof(dashmm_block_spawn_terminal_params_t) 
                        + payload_size;
    dashmm_block_spawn_terminal_params_t *input = malloc(input_size);
    assert(input != NULL);
    input->block_size = block_size;
    input->block_number = 0;
    input->payload_size = payload_size;
    memcpy(input->payload, payload, payload_size);
    
    hpx_call(base, terminal_action, input, input_size, retval);
    
    free(input);
  }
  
  return retval;
}


//
// Actions
//


int dashmm_block_spawn(void *args) {
  //interpret parameters
  dashmm_block_spawn_params_t *parms = args;
  
  //storage for the results
  hpx_addr_t results[2];
  results[0] = hpx_lco_future_new(parms->result_size);
  results[1] = hpx_lco_future_new(parms->result_size);
  assert(results[0] != HPX_NULL && results[1] != HPX_NULL);
  
  if (parms->last_block - parms->first_block == 1) {
    //get the left and right addresses
    hpx_addr_t targets[2];
    targets[0] = hpx_addr_add(parms->base, 
                                parms->block_size * parms->first_block,
                                parms->block_size);
    targets[1] = hpx_addr_add(parms->base, 
                                parms->block_size * parms->last_block,
                                parms->block_size);
    
    size_t input_size = sizeof(dashmm_block_spawn_terminal_params_t) 
                        + parms->payload_size;
    dashmm_block_spawn_terminal_params_t *inputs[2];
    inputs[0] = malloc(input_size);
    inputs[1] = malloc(input_size);
    assert(inputs[0] != NULL && inputs[1] != NULL);
    
    inputs[0]->block_size = parms->block_size;
    inputs[0]->block_number = parms->first_block;
    inputs[0]->payload_size = parms->payload_size;
    memcpy(inputs[0]->payload, parms->payload, parms->payload_size);
    
    inputs[1]->block_size = parms->block_size;
    inputs[1]->block_number = parms->last_block;
    inputs[1]->payload_size = parms->payload_size;
    memcpy(inputs[1]->payload, parms->payload, parms->payload_size);
    
    hpx_call(targets[0], parms->terminal_action, 
                        inputs[0], input_size, 
                        results[0]);
    hpx_call(targets[1], parms->terminal_action, 
                        inputs[1], input_size, 
                        results[1]);

    free(inputs[0]);
    free(inputs[1]);
  } else {
    size_t parms_size = sizeof(dashmm_block_spawn_params_t) 
                            + parms->payload_size;
    
    //we divide up the span
    dashmm_block_spawn_params_t *inputs[2];
    inputs[0] = malloc(parms_size);
    inputs[1] = malloc(parms_size);
    assert(inputs[0] != NULL && inputs[1] != NULL);
    
    memcpy(inputs[0], parms, parms_size);
    memcpy(inputs[1], parms, parms_size);
    
    int split = (parms->first_block + parms->last_block - 1) / 2;
    inputs[0]->last_block = split;
    inputs[1]->first_block = split + 1;
    //The division will put a single block chunk as the first
    assert(inputs[1]->last_block != inputs[1]->first_block);
    
    //get the left and right addresses
    hpx_addr_t targets[2];
    targets[0] = hpx_addr_add(parms->base, 
                                parms->block_size * parms->first_block,
                                parms->block_size);
    targets[1] = hpx_addr_add(parms->base, 
                                parms->block_size * (split + 1),
                                parms->block_size);
                               
    //launch the correct action
    if (inputs[0]->last_block == inputs[0]->first_block) {
      size_t tinput_size = sizeof(dashmm_block_spawn_terminal_params_t)
                            + parms->payload_size;
      dashmm_block_spawn_terminal_params_t *tinput;
      tinput = malloc(tinput_size);
      assert(tinput != NULL);
      
      tinput->block_size = parms->block_size;
      tinput->block_number = parms->first_block;
      tinput->payload_size = parms->payload_size;
      memcpy(tinput->payload, parms->payload, parms->payload_size);
      
      hpx_call(targets[0], parms->terminal_action, 
                    tinput, tinput_size, 
                    results[0]);

      free(tinput);
    } else {
      hpx_call(targets[0], dashmm_block_spawn_action, 
                inputs[0], parms_size, results[0]);
    }
                                
    hpx_call(targets[1], dashmm_block_spawn_action, 
                inputs[1], parms_size, results[1]);
    
    free(inputs[0]);
    free(inputs[1]);
  }
  
  
  //Wait for results, in whatever way is appropriate
  if (parms->result_size == 0) {
    //wait on results
    hpx_lco_wait_all(2, results, NULL);
      
    return HPX_SUCCESS;
  } else {
    //get the results
    void *bufs[2];
    bufs[0] = malloc(parms->result_size);
    bufs[1] = malloc(parms->result_size);
    assert(bufs[0] != NULL && bufs[1] != NULL);
    int szs[2] = {parms->result_size, parms->result_size};
     
    hpx_lco_get_all(2, results, szs, bufs, NULL);
      
    //reduce the results; after, bufs[0] will hold the reduced value
    parms->reduce(bufs[0], bufs[1]); 
    free(bufs[1]);
      
    //continue the value
    hpx_thread_continue_cleanup(parms->result_size, bufs[0], 
                                    dashmm_continue_cleanup, bufs[0]);
  }
}


int dashmm_point_volume(void *args) {
  hpx_addr_t points_gas = hpx_thread_current_target();
  dashmm_point_t *points;
  if (!hpx_gas_try_pin(points_gas, (void **)&points) ) {
    return HPX_RESEND;
  }
  
  //interpret arguments
  dashmm_block_spawn_terminal_params_t *parms = args;
  
  //how many points do we have?
  int n_points = parms->block_size / sizeof(dashmm_point_t);
  assert( (parms->block_size % sizeof(dashmm_point_t)) == 0);
  
  //set the starting volume
  dashmm_volume_t retval;
  retval.a[0] = points[0].pos[0];
  retval.a[1] = points[0].pos[1];
  retval.a[2] = points[0].pos[2];
  retval.b[0] = points[0].pos[0];
  retval.b[1] = points[0].pos[1];
  retval.b[2] = points[0].pos[2];
  
  //reduce over the points to get the volume
  for (int i = 1; i < n_points; ++i) {
    for (int dir = 0; dir < 3; ++dir) {
      if (points[i].pos[dir] < retval.a[dir]) {
        retval.a[dir] = points[i].pos[dir];
      }
      if (points[i].pos[dir] > retval.b[dir]) {
        retval.b[dir] = points[i].pos[dir];
      }
    }
  }
  
  hpx_gas_unpin(points_gas);
  
  //done
  hpx_thread_continue(sizeof(dashmm_volume_t), &retval);
}


int dashmm_tree_node_init(void *args) {
  hpx_addr_t node_gas = hpx_thread_current_target();
  dashmm_tree_node_t *node;
  if ( !hpx_gas_try_pin(node_gas, (void **)&node) ) {
    return HPX_RESEND;
  }
  
  //interpret parameters
  dashmm_block_spawn_terminal_params_t *wrap = args;
  dashmm_tree_node_init_params_t *parms = 
                    (dashmm_tree_node_init_params_t *)wrap->payload;
  
  //fill in the values
  int level = dashmm_tree_topnode_level(wrap->block_number);
  int index[3];
  dashmm_tree_topnode_index_in_level(wrap->block_number, level, index);
  node->vol = dashmm_tree_topnode_volume(parms->tree.vol, level, index);
  
#ifdef DEBUGID
  //make the ID
  node->ID = dashmm_tree_node_id(index, level);
#endif
  
  if (level) {
    int pindex[3];
    pindex[0] = index[0] / 2;
    pindex[1] = index[1] / 2;
    pindex[2] = index[2] / 2;
    node->parent = dashmm_tree_topnode_address(parms->tree.topnodes, 
                                                     level - 1, pindex);
  } else {
    node->parent = HPX_NULL;
  }
  
  node->n_points = 0;
  node->points = HPX_NULL;
  
  if (level < parms->tree.top_depth) {
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
    node->child[0] = hpx_lco_and_new(parms->point_block_count);
    assert(node->child[0] != HPX_NULL);
    
    node->child[1] = hpx_lco_future_new(0);
    assert(node->child[1] != HPX_NULL);
    
    node->child[2] = HPX_NULL;
    node->child[3] = HPX_NULL;
    node->child[4] = HPX_NULL;
    node->child[5] = HPX_NULL;
    node->child[6] = HPX_NULL;
    node->child[7] = HPX_NULL;
  }
  
  //done
  hpx_gas_unpin(node_gas);
  
  return HPX_SUCCESS;
}


int dashmm_tree_points_sort(void *args) {
  hpx_addr_t points_gas = hpx_thread_current_target();
  dashmm_point_t *points;
  if ( !hpx_gas_try_pin(points_gas, (void **)&points) ) {
    return HPX_RESEND;
  }

  //get parameters
  dashmm_block_spawn_terminal_params_t *wrap = args;
  dashmm_tree_t *tree = (dashmm_tree_t *)wrap->payload;
  
  //how many points?
  int n_points = wrap->block_size / sizeof(dashmm_point_t);
  
  //set up the bins for counting
  int n_bins = _onlevels[tree->top_depth];
  //int *bins = calloc(sizeof(int), n_bins);
  int *bins = malloc(sizeof(int) * n_bins);
  assert(bins != NULL);
  for (int ibin = 0; ibin < n_bins; ++ibin) {
    bins[ibin] = 0;
  }

  //count for each topnode; mark bins for the points
  for (int ipt = 0; ipt < n_points; ++ipt) {
    //figure out my on level index
    int index[3];
    dashmm_tree_topnode_from_point(tree->vol, tree->top_depth, 
                                    &points[ipt], index);
    //figure out linear offset from that
    int mybin = dashmm_tree_topnode_offset_on_level(tree->top_depth, index);
    assert(mybin < n_bins);
    //bin it
    bins[mybin] += 1;
    points[ipt].sort = mybin;
  }
   
  //tell topnodes the results
  hpx_addr_t *results = malloc(sizeof(hpx_addr_t) * n_bins);
  assert(results != NULL);
  for (int ires = 0; ires < n_bins; ++ires) {
    results[ires] = hpx_lco_future_new(sizeof(hpx_addr_t));
    assert(results[ires] != HPX_NULL);
  }
  for (int inode = 0; inode < n_bins; ++inode) {
    int total = _offsets[tree->top_depth] + inode;
    hpx_addr_t target = hpx_addr_add(tree->topnodes, 
                                     total * sizeof(dashmm_tree_node_t), 
                                     sizeof(dashmm_tree_node_t));
    hpx_call(target, dashmm_tree_inform_topnode_action, 
             &bins[inode], sizeof(int), results[inode]);
  }
  
  //perform the sort
  dashmm_point_bin_sort(points, n_points, bins, n_bins);

  //wait on return from topnode action
  hpx_lco_wait_all(n_bins, results, NULL);
  
  //hpx_gas_memput to move points into correct location
  hpx_addr_t sendsdone = hpx_lco_and_new(n_bins);
  int offset = 0;
  for (int ibin = 0; ibin < n_bins; ++ibin) {
    //Only do this if there are points in the given bin
    if (bins[ibin]) {
      //get address of destination
      hpx_addr_t dest;
      hpx_lco_get(results[ibin], sizeof(dest), &dest);
      //compute amount to send
      size_t cpysize = sizeof(dashmm_point_t) * bins[ibin];
      //send it
      hpx_gas_memput(dest, &points[offset], cpysize, HPX_NULL, sendsdone);
      //update offset
      offset += bins[ibin];
    } else {
      //otherwise, just set the send counter
      hpx_lco_and_set(sendsdone, HPX_NULL);
    }
  }
  hpx_lco_wait(sendsdone);
  hpx_lco_delete(sendsdone, HPX_NULL);
  
  //done
  for (int ibin = 0; ibin < n_bins; ++ibin) {
    hpx_lco_delete(results[ibin], HPX_NULL);
  }
  free(bins);
  free(results);
  hpx_gas_unpin(points_gas);
  
  return HPX_SUCCESS;
}


int dashmm_tree_inform_topnode(void *args) {
  hpx_addr_t node_gas = hpx_thread_current_target();
  dashmm_tree_node_t *node;
  if ( !hpx_gas_try_pin(node_gas, (void **)&node) ) {
    return HPX_RESEND;
  }
  
  //interpret the parameter
  int in_count = *((int *)args);
  
  //atomically increment the counter
  int myoffset = sync_fadd(&(node->n_points), in_count, SYNC_ACQ_REL);
  
  //set the and gate
  hpx_lco_and_set(node->child[0], HPX_NULL);
  
  //wait on the allocation signal
  hpx_lco_wait(node->child[1]);
  
  //compute the return value
  hpx_addr_t result = hpx_addr_add(node->points, 
                                   sizeof(dashmm_point_t) * myoffset,
                                   sizeof(dashmm_point_t) * node->n_points);
  
  //done
  hpx_gas_unpin(node_gas);
  
  hpx_thread_continue(sizeof(hpx_addr_t), &result);
}


int dashmm_tree_node_alloc(void *args) {
  hpx_addr_t node_gas = hpx_thread_current_target();
  dashmm_tree_node_t *node;
  if ( !hpx_gas_try_pin(node_gas, (void **)&node) ) {
    return HPX_RESEND;
  }
  
  //get parameters
  dashmm_block_spawn_terminal_params_t *wrap = args;
  dashmm_tree_node_alloc_params_t *parms = 
                    (dashmm_tree_node_alloc_params_t *)wrap->payload;
  
  //wait for the incoming counts
  hpx_lco_wait(node->child[0]);

  
  //allocate space for incoming points
  if (node->n_points) {
    node->points = hpx_gas_alloc(sizeof(dashmm_point_t) * node->n_points);
    assert(node->points != HPX_NULL);
  }
  
  //set allocated signal
  hpx_lco_set(node->child[1], 0, NULL, HPX_NULL, HPX_NULL);
  
  //wait for sorting to finish
  hpx_lco_wait(parms->sortdone);
  
  //clean up some stuff
  hpx_lco_delete(node->child[0], HPX_NULL);
  node->child[0] = HPX_NULL;
  hpx_lco_delete(node->child[1], HPX_NULL);
  node->child[1] = HPX_NULL;
  
  //begin refinement of this topnode
  if (node->n_points) {
    dashmm_tree_node_refine_params_t input;
    input.refine_limit = parms->refine_limit;
    hpx_call_sync(node_gas, dashmm_tree_node_refine_action, 
                &input, sizeof(input), 
                NULL, 0);
  }
  
  //done
  hpx_gas_unpin(node_gas);
  
  return HPX_SUCCESS;
}


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
#ifdef DEBUGID
        input.ID = dashmm_tree_node_index_for_child(node->ID, i);
#endif
        
        hpx_call(node->child[i], dashmm_tree_init_child_action,
                 &input, sizeof(input), alldone);
      }
    }
    
    hpx_lco_wait(alldone);
    hpx_lco_delete(alldone, HPX_NULL);
    
  } 
  //else {
    //we are happy how we are
    //FOR NOW: nothing here, but in the future here is where we can start
    // the moment computation
    //FUTURE: We might want to create a local copy of the point data, so that
    // some operations can proceed quickly. Or not.
    
    //THIS branch is also executed when the node has no points. This is
    // relevant for topnodes that did not end up with any points. In that case,
    // we do not want this to do anything.
  //}
  
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
  
  //Sort the points
  dashmm_point_bin_sort(points, parms->n_points, result.n_points, 8);
  
  //setup offsets for the different blocks
  int total_offset = 0;
  for(int i = 0; i < 8; ++i) {
    result.offsets[i] = hpx_addr_add(points_gas, 
                                    total_offset * sizeof(dashmm_point_t),
                                    parms->n_points * sizeof(dashmm_point_t));
    total_offset += result.n_points[i];
  }
  assert(total_offset == parms->n_points);
  
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
#ifdef DEBUGID
  node->ID = parms->ID;
#endif
  
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


