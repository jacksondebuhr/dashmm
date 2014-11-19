// =============================================================================
//  DASHMM 
//
//  Test driver for draft of tree construction for DASHMM
//  test-driver.c
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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "hpx/hpx.h"

#include "dashmm-tree.h"


//define this to test the branch only construction - in SMP
//#define BRANCHESONLY


hpx_action_t hpx_main_action;
hpx_action_t create_point_data_action;

typedef struct {
  int n_points;
  int refinement_limit;
  int top_depth;
} hpx_main_params_t;

typedef struct {
  int n_points;
  int n_per_block;
  int n_blocks;
  int i_block;
  //eventually some option for data type
  // (i.e. cube; sphere, two cubes, two spheres)
} create_point_data_params_t;

int hpx_main(void *args);
int create_point_data(void *args);


static void usage(FILE *f) {
  fprintf(f, "Usage: ygg2 [options] <number of points> <refinement limit> <top depth> \n"
          "\t-c, number of cores to run on\n"
          "\t-t, number of scheduler threads\n"
          "\t-T, select a transport by number (see hpx_config.h)\n"
          "\t-D, all localities wait for debugger\n"
          "\t-d, wait for debugger at specific locality\n"
          "\t-l, set logging level\n"
          "\t-s, set stack size\n"
          "\t-p, set per-PE global heap size\n"
          "\t-r, set send/receive request limit\n"
          "\t-h, show help\n");
}


int main(int argc, char **argv) {
  //initialize the runtime
  if( hpx_init(&argc, &argv) ) {
    fprintf(stderr, "HPX failed to initialize.\n");
    return 1;
  }

  int opt = 0;
  while ((opt = getopt(argc, argv, "c:t:T:d:Dl:s:p:r:h")) != -1) {
    switch (opt) {
    case 'c':
      break;
    case 't':
      break;
    case 'T':
      break;
    case 'D':
      break;
    case 'd':
      break;
    case 'l':
      break;
    case 's':
      break;
    case 'p':
      break;
    case 'r':
      break;
    case 'h':
      usage(stdout);
      return 0;
    case '?':
    default:
      usage(stderr);
      return -1;
    }
  }
  
  argc -= optind;
  argv += optind;

  if(argc < 3) {
    usage(stderr);
    return -1;
  }
  
  int n_points = atoi(argv[argc-3]);
  int refinement_limit = atoi(argv[argc-2]);
  int top_depth = atoi(argv[argc-1]);
  
  
  
  //register actions
  hpx_main_action = HPX_REGISTER_ACTION(hpx_main);
  create_point_data_action = HPX_REGISTER_ACTION(create_point_data);
  dashmm_tree_register_actions();
  
  //setup parameters to main action
  hpx_main_params_t input;
  input.n_points = n_points;
  input.refinement_limit = refinement_limit;
  input.top_depth = top_depth;
  
  //start it up
  return hpx_run(hpx_main_action, &input, sizeof(input));
}


int hpx_main(void *args) {
  //interpret params
  hpx_main_params_t *parms = args;
  
  
  //run an action to create some point data
  create_point_data_params_t input;
  input.n_points = parms->n_points;
  //more options eventually...
  
#ifdef BRANCHESONLY
  input.n_per_block = parms->n_points;
  input.n_blocks = 1;
  input.i_block = 0;
  hpx_addr_t points = hpx_gas_alloc(sizeof(dashmm_point_t) * input.n_points);
#else
  int numlocs = hpx_get_num_ranks();
  int numperblock = input.n_points / numlocs;
  if(input.n_points % numlocs != 0) {
    ++numperblock;
  }
  input.n_points = numlocs * numperblock;
  input.n_per_block = numperblock;
  input.n_blocks = numlocs;
  input.i_block = 0;
  hpx_addr_t points = hpx_gas_global_alloc(numlocs, 
                                    sizeof(dashmm_point_t) * numperblock);
#endif
  assert(points != HPX_NULL);

fprintf(stdout, "Beginning point creation\n");fflush(stdout);  

  hpx_call_sync(points, create_point_data_action, 
                &input, sizeof(input),
                NULL, 0);

fprintf(stdout, "Done with point creation\n");fflush(stdout);
  
  
  //start a timer
  hpx_time_t start_time = hpx_time_now();
  
  dashmm_volume_t rootvol;
  rootvol.a[0] = 0.0;
  rootvol.a[1] = 0.0;
  rootvol.a[2] = 0.0;
  rootvol.b[0] = 1.0;
  rootvol.b[1] = 1.0;
  rootvol.b[2] = 1.0;
  
#ifdef BRANCHESONLY  
  //create the tree
  //NOTE that is only using a subset of the eventual functionality
  hpx_addr_t root = hpx_gas_alloc(sizeof(dashmm_tree_node_t));
  assert(root != HPX_NULL);
  
  dashmm_invoke_tree_init_child_sync(root, rootvol, HPX_NULL, 
                points, parms->n_points, parms->refinement_limit);
#else
  //Do some stuff for a full tree...
  hpx_addr_t tree = dashmm_tree_create(points, parms->n_points, 
                        input.n_per_block, rootvol, 
                        parms->refinement_limit,
                        parms->top_depth);
#endif  
  
  //stop the timer
  hpx_time_t end_time = hpx_time_now();
  double dt_ms = hpx_time_diff_ms(start_time, end_time);
  
  
  //put out some information
  fprintf(stdout, "%d %lg\n", parms->n_points, dt_ms);
  //more stuff that will be interesting I guess would be the number of nodes
  // the number of cores used and so on
  
  hpx_shutdown(HPX_SUCCESS);
}


int create_point_data(void *args) {
  //Please NOTE: This is a  more or less silly way to do the synchronization
  // for this. However, this is just the initialization, so who cares.

  hpx_addr_t points_gas = hpx_thread_current_target();
  dashmm_point_t *points;
  if( !hpx_gas_try_pin(points_gas, (void **)&points) ) {
    return HPX_RESEND;
  }

  //read in the options
  create_point_data_params_t *parms = args;
  
  //put out message about number of points and data set type
  //fprintf(stdout, 
  //        "Creating point data for %d points distributed in a unit cube\n",
  //        parms->n_points);
  
  //start a basic RNG
  srand(time(NULL));
  
  //fill the data
  for(int i = 0; i < parms->n_per_block; ++i) {
    points[i].pos[0] = ((double)rand())/((double)RAND_MAX);
    points[i].pos[1] = ((double)rand())/((double)RAND_MAX);
    points[i].pos[2] = ((double)rand())/((double)RAND_MAX);
  }
  
  hpx_gas_unpin(points_gas);
  
  
  //decide if there is a next block to handle
  parms->i_block += 1;
  if(parms->i_block < parms->n_blocks) {
    //send the action to the next block
    hpx_addr_t target = hpx_addr_add(points_gas, 
                            sizeof(dashmm_point_t) * parms->n_per_block,
                            sizeof(dashmm_point_t) * parms->n_per_block);
    hpx_call_sync(target, create_point_data_action, 
                    parms, sizeof(*parms), NULL, 0);
  }
  
  
  return HPX_SUCCESS;
}


