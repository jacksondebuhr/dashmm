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

#include "particle.h"
#include "tree.h"


void print_usage(FILE *fd, char *prog) {
  fprintf(fd, "Usage: %s <N parts> <Partition Limit> <theta> <domain size>\n",
          prog);
}


hpx_addr_t generate_parts(int n_parts, double domain_length, hpx_addr_t where) {
  //generate n_parts over the given domain length
  hpx_addr_t parts_gas = hpx_gas_alloc_local_at_sync(1,
                            sizeof(Particle) * n_parts, 0, where);
  assert(parts_gas != HPX_NULL);
  Particle *parts = NULL;
  assert(hpx_gas_try_pin(parts_gas, (void **)&parts));

  for (int i = 0; i < n_parts; ++i) {
    parts[i].pos = ((double)rand() / (double)RAND_MAX) * domain_length;
    parts[i].mass = ((double)rand() / (double)RAND_MAX) * 0.9 + 0.1;
  }

  hpx_gas_unpin(parts_gas);

  return parts_gas;
}


int hpx_main(int n_parts, int n_partition, double theta_c,
             double domain_size) {
  broadcast_domain_size(domain_size);

  hpx_addr_t root = create_node(0.0, domain_size);
  hpx_addr_t parts = generate_parts(n_parts, domain_size, root);

  hpx_time_t t_start = hpx_time_now();

  partition_node_sync(root, parts, n_parts, n_partition);

  hpx_addr_t alldone = hpx_lco_and_new(n_parts);
  assert(alldone != HPX_NULL);
  spawn_potential_computation(root, alldone, theta_c);
  hpx_lco_wait(alldone);
  hpx_lco_delete_sync(alldone);

  hpx_time_t t_end = hpx_time_now();

  fprintf(stdout, "Time to compute %d potentials across %d localities:\n",
            n_parts, hpx_get_num_ranks());
  fprintf(stdout, "      %lg [ms]\n", hpx_time_diff_ms(t_start, t_end));

  hpx_exit(0, NULL);
}
HPX_ACTION(HPX_DEFAULT, 0, hpx_main_action, hpx_main,
           HPX_INT, HPX_INT, HPX_DOUBLE, HPX_DOUBLE);


int main(int argc, char **argv) {
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
