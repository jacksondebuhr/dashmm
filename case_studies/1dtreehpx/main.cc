#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <vector>
#include <algorithm>
#include <random>
#include <functional>

#include "hpx/hpx.h"

#include "tree.h"
#include "particle.h"


void print_usage(FILE *fd, char *prog) {
  fprintf(fd, "Usage: %s <N parts> <Partition limit>"
              " <theta limit> <domain size>\n", prog);
}


hpx_addr_t generate_parts(int n_parts, double domain_length) {
  hpx_addr_t parts_gas = hpx_gas_alloc_local_at_sync(
                            sizeof(Particle) * n_parts, 0, HPX_HERE);
  assert(parts_gas != HPX_NULL);
  Particle *parts{nullptr};
  //may need to do static_cast<void **>(&parts) instead
  assert(hpx_gas_try_pin(parts_gas, (void **)&parts));

  std::mt19937 engine;
  std::uniform_real_distribution<double> posdistrib(0.0, domain_length);
  std::uniform_real_distribution<double> massdistrib(0.1, 1.0);
  auto xroll = std::bind(posdistrib, engine);
  auto mroll = std::bind(massdistrib, engine);

  for (int i = 0; i < n_parts; ++i) {
    double x = xroll();
    double m = mroll();
    parts[i] = Particle(x, m);
  }
  
  std::sort(&parts[0], &parts[n_parts]);
  for (int i = 0; i < n_parts; ++i) {
    hpx_addr_t data = hpx_addr_add(parts_gas, sizeof(Particle) * i, 
                                   sizeof(Particle) * n_parts);
    parts[i].set_data(data);
  }
  
  hpx_gas_unpin(parts_gas);
  
  return parts_gas;
}


int hpx_main(int n_parts, int n_partition, 
             double theta_c, double domain_size) {
  //generate the particles
  hpx_addr_t parts = generate_parts(n_parts, domain_size);
  Particle *P {nullptr};
  assert(hpx_gas_try_pin(parts, (void **)&P));
  
  //find the overall domain size
  double min{P[0].x()};
  double max{P[n_parts - 1].x()};
  double center{0.5 * (min + max)};
  double span{0.5001 * (max - min)};
  double upper{center + span};
  double lower{center - span};
  
  //construct the tree
  Node root(lower, upper);
  root.partition_sync(parts, 0, n_parts, n_partition);
  
  //compute moments
  root.compute_moment_sync();
  
  //loop over particles and compute potential for each
  hpx_time_t t1 = hpx_time_now();
  hpx_addr_t alldone = hpx_lco_and_new(n_parts);
  for (int i = 0; i < n_parts; ++i) {
    double pos{P[i].x()};
    root.compute_potential_and_continue(pos, i, theta_c, 
                                        P[i].set_approx_cont(),
                                        alldone);
  }
  hpx_lco_wait(alldone);
  hpx_time_t t2 = hpx_time_now();
  hpx_lco_delete(alldone, HPX_NULL);
  double dt = hpx_time_diff_ms(t1, t2);
  fprintf(stdout, "It took %lg ms to compute %d potentials\n",
                    dt, n_parts);
  
  //done
  //NOTE: We explicitly call the destructor, since HPX-5 does not actually
  // unwind the stack, and destroy stack objects.
  hpx_gas_unpin(parts);
  root.~Node();
  hpx_gas_free(parts, HPX_NULL);
  
  hpx_shutdown(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, 0, 
           hpx_main_action, hpx_main,
           HPX_INT, HPX_INT, HPX_DOUBLE, HPX_DOUBLE);


int main(int argc, char **argv) {
  if (hpx_init(&argc, &argv)) {
    hpx_print_help();
    return -1;
  }
  
  //Interpret the command line arguments
  if (argc != 5) {
    print_usage(stderr, argv[0]);
    return 0;
  }
  
  int n_parts{atoi(argv[1])};
  int n_partition{atoi(argv[2])};
  double theta_c{atof(argv[3])};
  double domain_size{atof(argv[4])};
  
  return hpx_run(&hpx_main_action, &n_parts, &n_partition, 
                 &theta_c, &domain_size);
}


