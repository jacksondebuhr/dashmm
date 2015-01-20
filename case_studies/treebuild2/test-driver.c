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
#include <math.h>

#include "hpx/hpx.h"

#include "dashmm-tree.h"



hpx_action_t hpx_main_action;
hpx_action_t create_point_data_action;

typedef struct {
  int set_type;
  int n_points;
  int refinement_limit;
  int top_depth;
  int block_count;
} hpx_main_params_t;

typedef struct {
  int set_type;
  int n_points;
  int n_per_block;
  int n_blocks;
  int i_block;
} create_point_data_params_t;

int hpx_main(void *args);
void create_point_data(hpx_addr_t points_gas, int n_per_block, int n_blocks, 
                        int set_type);


int main(int argc, char **argv) {
  if (argc < 6) {
    fprintf(stdout, "Usage: %s <point set type> <n points> <refinement limit> <top depth> <number of point blocks>\n", argv[0]);
    fprintf(stdout, 
     "Valid inputs will have n_points divisible by the number of point blocks\n");
    fprintf(stdout, "<point set type>\n");
    fprintf(stdout, "\t0 - uniform cube\n\t1 - Hernquist Halo\n");
    fprintf(stdout, "\t2 - disjoint spheres\n");
    return 0;
  }
  
  int set_type = atoi(argv[1]);
  int n_points = atoi(argv[2]);
  int refinement_limit = atoi(argv[3]);
  int top_depth = atoi(argv[4]);
  int block_count = atoi(argv[5]);
  
  
  //initialize the runtime
  if ( hpx_init(&argc, &argv) ) {
    fprintf(stderr, "HPX failed to initialize.\n");
    return 1;
  }
  
  
  //register actions
  HPX_REGISTER_ACTION(hpx_main, &hpx_main_action);
  dashmm_tree_register_actions();
  
  //setup parameters to main action
  hpx_main_params_t input;
  input.set_type = set_type;
  input.n_points = n_points;
  input.refinement_limit = refinement_limit;
  input.top_depth = top_depth;
  input.block_count = block_count;
  
  //start it up
  return hpx_run(&hpx_main_action, &input, sizeof(input));
}


int hpx_main(void *args) {
  //interpret params
  hpx_main_params_t *parms = args;

  //compute some numbers
  int numblocks = parms->block_count * hpx_get_num_ranks();
  int numperblock = parms->n_points / numblocks;
  if (parms->n_points % numblocks != 0) {
    ++numperblock;
  }

  if (parms->n_points % numblocks) {
    fprintf(stderr, 
      "Number of points must be evenly divisible by the number of blocks.\n");
    return HPX_SUCCESS;
  }

//fprintf(stdout, "Trying to get numblocks=%d with %d parts per block of %lu bytes each\n", numblocks, numperblock, sizeof(dashmm_point_t)*numperblock);


  hpx_addr_t points = hpx_gas_global_alloc(numblocks, 
                                    sizeof(dashmm_point_t) * numperblock);

  //assert(points != HPX_NULL);

  create_point_data(points, numperblock, numblocks, parms->set_type);
  
  //start a timer
  hpx_time_t start_time = hpx_time_now();
  
  //Do some stuff for a full tree...
  hpx_addr_t tree = dashmm_tree_create(points, parms->n_points, 
                        numperblock, 
                        parms->refinement_limit,
                        parms->top_depth);
  
  //stop the timer
  hpx_time_t end_time = hpx_time_now();
  double dt_ms = hpx_time_diff_ms(start_time, end_time);
  
  //put out some information
  fprintf(stdout, "%d %lg\n", parms->n_points, dt_ms);
  //more stuff that will be interesting I guess would be the number of nodes
  // the number of cores used and so on
  
  //free up the tree
  dashmm_tree_destroy(tree);
  
  hpx_shutdown(HPX_SUCCESS);
}


void create_point_data(hpx_addr_t points_gas, int n_per_block, int n_blocks, 
                            int set_type) {
  //start a basic RNG
  srand(time(NULL));

  //create local space for the points
  dashmm_point_t *points = malloc(sizeof(dashmm_point_t) * n_per_block);
  assert(points != NULL);
  
  for (int iblock = 0; iblock < n_blocks; ++iblock) {
    //fill the data
    switch (set_type) {
    case 0:
      //Unit cube
      for (int i = 0; i < n_per_block; ++i) {
        points[i].pos[0] = ((double)rand())/((double)RAND_MAX);
        points[i].pos[1] = ((double)rand())/((double)RAND_MAX);
        points[i].pos[2] = ((double)rand())/((double)RAND_MAX);
      }
      break;
    case 1:
      for (int i = 0; i < n_per_block; ++i) {
        //Hernquist Halo
        double x = 25.0 / 36.0 * ((double)rand())/((double)RAND_MAX);
        double theta = 3.1415926535 * ((double)rand())/((double)RAND_MAX);
        double phi = 2.0 * 3.1415926535 * ((double)rand())/((double)RAND_MAX);
        double r = 0.1 * x / (1.0 - x) * (1.0 + sqrt(1.0 / x)) * 20.0;
        points[i].pos[0] = 0.5 + r * sin(theta) * cos(phi);
        points[i].pos[1] = 0.5 + r * sin(theta) * sin(phi);
        points[i].pos[2] = 0.5 + r * cos(theta);
      }
      break;
    case 2:
      for (int i = 0; i < n_per_block; ++i) {
        //Two disjoint spheres
        double r = ((double)rand())/((double)RAND_MAX);
        double theta = 3.1415926535 * ((double)rand())/((double)RAND_MAX);
        double phi = 2.0 * 3.1415926535 * ((double)rand())/((double)RAND_MAX);
        double xcen[3] = {0.0};
        if (rand() % 2) {        //point for sphere 1
          r *= 0.5;
          xcen[0] = 0.5;
          xcen[1] = 0.5;
          xcen[2] = 0.5;
        } else {                //point for sphere 2
          r *= 0.5;
          xcen[0] = -0.5;
          xcen[1] = -0.5;
          xcen[2] = -0.5;
        }
        points[i].pos[0] = xcen[0] + r * sin(theta) * cos(phi);
        points[i].pos[1] = xcen[1] + r * sin(theta) * sin(phi);
        points[i].pos[2] = xcen[2] + r * cos(theta);
      }
      break;
    default:
      //Unit cube
      for (int i = 0; i < n_per_block; ++i) {
        points[i].pos[0] = ((double)rand())/((double)RAND_MAX);
        points[i].pos[1] = ((double)rand())/((double)RAND_MAX);
        points[i].pos[2] = ((double)rand())/((double)RAND_MAX);
      }
      break;
    }
    
    //now do a memput
    hpx_addr_t target = hpx_addr_add(points_gas, 
                           sizeof(dashmm_point_t) * n_per_block * iblock,
                           sizeof(dashmm_point_t) * n_per_block);

    hpx_addr_t local = hpx_lco_future_new(0);
    hpx_addr_t global = hpx_lco_future_new(0);
    assert(local != HPX_NULL && global != HPX_NULL);
    
    hpx_gas_memput(target, points, sizeof(dashmm_point_t) * n_per_block,
                        local, global);
                        
    hpx_lco_wait(local);
    hpx_lco_wait(global);
    hpx_lco_delete(local, HPX_NULL);
    hpx_lco_delete(global, HPX_NULL);
  }

}

