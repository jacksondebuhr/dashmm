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

#include <cstdio>
#include <cstdlib>

#include "hpx/hpx.h"

#include "tree.h"


void print_usage(FILE *fd, char *prog) {
  fprintf(fd, "Usage: %s <N parts> <Partition Limit> <theta> <domain size>\n",
          prog);
}


Particle *generate_parts(int n_parts, double domain_length) {
  Particle *parts = new Particle[n_parts];

  for (int i = 0; i < n_parts; ++i) {
    parts[i].pos = ((double)rand() / (double)RAND_MAX) * domain_length;
    parts[i].mass = ((double)rand() / (double)RAND_MAX) * 0.9 + 0.1;
  }

  return parts;
}


int hpx_main(int n_parts, int n_partition, double theta_c,
             double domain_size) {
  Node *root = new Node(0.0, domain_size);
  Particle *parts = generate_parts(n_parts, domain_size);


  // Create the tree
  root->partition(parts, n_parts, n_partition);
  fprintf(stdout, "Done with partitioning!\n");


  root->compute_moments();
  fprintf(stdout, "Done computing moments.\n");


  // Destroy the tree
  delete root;
  fprintf(stdout, "Done with tree deletion.\n");


  // Destroy particles
  delete [] parts;

  hpx_exit(0, NULL);
}
HPX_ACTION(HPX_DEFAULT, 0, hpx_main_action, hpx_main,
           HPX_INT, HPX_INT, HPX_DOUBLE, HPX_DOUBLE);


int main(int argc, char **argv) {
  Node::register_actions();

  if (hpx_init(&argc, &argv)) {
    hpx_print_help();
    return -1;
  }

  if (argc != 5) {
    print_usage(stderr, argv[0]);
    return 0;
  }

  int n_parts = atoi(argv[1]);
  int n_partition = atoi(argv[2]);
  double theta_c = atof(argv[3]);
  double domain_size = atof(argv[4]);


  int err = hpx_run(&hpx_main_action, NULL, &n_parts, &n_partition,
                    &theta_c, &domain_size);

  hpx_finalize();

  return err;
}
