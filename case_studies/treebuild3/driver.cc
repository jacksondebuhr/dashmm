#include <iostream>
#include <getopt.h>

#include "dashmm/array.h"

#include "tree.h"


void usage(std::string exec) {
  std::cout << "Usage: " << exec << " --scaling=[w/s] "
            << "--datatype=[c/s] "
            << "--nsrc=num " << "--ntar=num "
            << "--threshold=num " << "--nseed=num " << std::endl;
  std::cout << "Note: --nseed=num required if --scaling=s" << std::endl;
}

/// The main action of the code starts here
int main_handler(char scaling, char datatype,
                 int nsrc, int ntar, int threshold, int nseed,
                 hpx_addr_t source_gas, hpx_addr_t target_gas) {
  int num_ranks = hpx_get_num_ranks();

  // Partition points to create dual tree

  hpx_time_t timer_start = hpx_time_now();

  // Determine domain geometry
  hpx_addr_t domain_geometry =
    hpx_lco_reduce_new(num_ranks, sizeof(double) * 6,
                       domain_geometry_init_action,
                       domain_geometry_op_action);

  hpx_bcast_lsync(set_domain_geometry_action, HPX_NULL,
                  &source_gas, &target_gas, &domain_geometry);

  double var[6];
  hpx_lco_get(domain_geometry, sizeof(double) * 6, &var);
  hpx_lco_delete_sync(domain_geometry);
  double size_x = var[1] - var[0];
  double size_y = var[3] - var[2];
  double size_z = var[5] - var[4];
  size = fmax(size_x, fmax(size_y, size_z));
  corner_x = (var[1] + var[0] - size) / 2;
  corner_y = (var[3] + var[2] - size) / 2;
  corner_z = (var[5] + var[4] - size) / 2;

  // Measure the reduction time
  hpx_time_t timer_domain = hpx_time_now();

  // Choose uniform partition level such that the number of grids is no less
  // than the number of ranks
  unif_level = ceil(log(num_ranks) / log(8)) + 1;

  int dim3 = pow(8, unif_level);
  unif_count = hpx_lco_reduce_new(num_ranks, sizeof(int) * (dim3 * 2),
                                  int_sum_ident_op,
                                  int_sum_op);

  hpx_bcast_rsync(init_partition_action, &unif_count, &unif_level, &threshold,
                  &corner_x, &corner_y, &corner_z, &size);

  // Measure the partitioning setup. This is all just getting this and that
  // ready to do. No real work is done in the above.
  hpx_time_t timer_middle = hpx_time_now();

  hpx_bcast_rsync(create_dual_tree_action, &source_gas, &target_gas);

  // All done, spit out some timing information before cleaning up and halting
  hpx_time_t timer_end = hpx_time_now();
  double elapsed_total = hpx_time_diff_ms(timer_start, timer_end) / 1e3;
  double elapsed_reduction = hpx_time_diff_ms(timer_start, timer_domain) / 1e3;
  double elapsed_first = hpx_time_diff_ms(timer_domain, timer_middle) / 1e3;
  double elapsed_second = hpx_time_diff_ms(timer_middle, timer_end) / 1e3;
  std::cout << "Dual tree creation time: " << elapsed_total << "\n";
  std::cout << "  Reduction: " << elapsed_reduction << "\n";
  std::cout << "  Setup: " << elapsed_first << "\n";
  std::cout << "  Partition: " << elapsed_second << "\n";

  // This is just deleting the allocated resources. Nothing too special here.
  hpx_bcast_rsync(finalize_partition_action, NULL, 0);

  hpx_lco_delete_sync(unif_count);

  hpx_exit(0, nullptr);
}
HPX_ACTION(HPX_DEFAULT, 0, main_action, main_handler,
           HPX_CHAR, HPX_CHAR, HPX_INT, HPX_INT, HPX_INT, HPX_INT,
           HPX_ADDR, HPX_ADDR);

int main(int argc, char *argv[]) {
  // default values
  char scaling = 'w'; // strong scale, 'w' for weak scaling
  char datatype = 'c'; // cubic point distribution, 's' for sphere distrbution
  int nsrc = 20; // total number of sources for strong scaling
                 // number of sources per rank for weak scaling
  int ntar = 20;
  int threshold = 1;
  int nseed = -1; // Number of seeds used for generating input for strong
                  // scaling test

  if (hpx_init(&argc, &argv)) {
    std::cout << "HPX: failed to initialize" << std::endl;
  }

  // Parse command line
  int opt = 0;
  static struct option long_options[] = {
    {"scaling", required_argument, 0, 's'},
    {"datatype", required_argument, 0, 'd'},
    {"nsrc", required_argument, 0, 'n'},
    {"ntar", required_argument, 0, 'm'},
    {"threshold", required_argument, 0, 'l'},
    {"nseed", required_argument, 0, 'r'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
  };

  int long_index = 0;
  bool valid_arguments = true;
  while ((opt = getopt_long(argc, argv, "s:d:e:n:m:l:r:h",
                            long_options, &long_index)) != -1) {
    switch (opt) {
    case 's':
      scaling = *optarg;
      break;
    case 'd':
      datatype = *optarg;
      break;
    case 'n':
      nsrc = atoi(optarg);
      break;
    case 'm':
      ntar = atoi(optarg);
      break;
    case 'l':
      threshold = atoi(optarg);
      break;
    case 'r':
      nseed = atoi(optarg);
      break;
    case 'h':
      usage(std::string{argv[0]});
      valid_arguments = false;
      break;
    }
  }

  if (scaling == 's') {
    // For strong scaling test, @p nseed must be overwritten.
    if (nseed <= 0) {
      valid_arguments = false;
      usage(std::string{argv[0]});
    }
  }

  if (!valid_arguments) {
    hpx_finalize();
    return -1;
  }

  // First get the data setup
  dashmm::Array<Point> sources = generate_points(scaling, datatype,
                                                  nsrc, nseed, 0);
  dashmm::Array<Point> targets = generate_points(scaling, datatype,
                                                  ntar, nseed, nseed);

  // Then build the tree
  // TODO - make this actually work
  hpx_addr_t ssend = sources.data();  // TODO: this will eventually disappear
  hpx_addr_t tsend = targets.data();  // behind an interface
  if (hpx_run(&main_action, nullptr, &scaling, &datatype,
              &nsrc, &ntar, &threshold, &nseed, &ssend, &tsend)) {
    std::cout << "Failed to run main action" << std::endl;
  }

  // Then delete the tree

  // Then delete the data
  sources.destroy();
  targets.destroy();

  hpx_finalize();

  return 0;
}
