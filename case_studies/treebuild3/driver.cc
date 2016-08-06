#include <iostream>
#include <getopt.h>

#include "dashmm/array.h"
using dashmm::Array;

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
  Array<Point> sources{source_gas};
  Array<Point> targets{target_gas};

  // Partition points to create dual tree

  hpx_time_t timer_start = hpx_time_now();

  // Determine domain geometry
  compute_domain_geometry(sources, targets);
  hpx_time_t timer_domain = hpx_time_now();

  // Perform some basic setup
  setup_basic_data(threshold);
  hpx_time_t timer_middle = hpx_time_now();

  // Start the partitioning proper
  create_dual_tree(sources, targets);
  hpx_time_t timer_end = hpx_time_now();

  // All done, spit out some timing information before cleaning up and halting
  double elapsed_total = hpx_time_diff_ms(timer_start, timer_end) / 1e3;
  double elapsed_reduction = hpx_time_diff_ms(timer_start, timer_domain) / 1e3;
  double elapsed_first = hpx_time_diff_ms(timer_domain, timer_middle) / 1e3;
  double elapsed_second = hpx_time_diff_ms(timer_middle, timer_end) / 1e3;
  std::cout << "Dual tree creation time: " << elapsed_total << "\n";
  std::cout << "  Reduction: " << elapsed_reduction << "\n";
  std::cout << "  Setup: " << elapsed_first << "\n";
  std::cout << "  Partition: " << elapsed_second << "\n";

  // Here we clean up the allocated resources
  finalize_partition();

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
  Array<Point> sources = generate_points(scaling, datatype,
                                                  nsrc, nseed, 0);
  Array<Point> targets = generate_points(scaling, datatype,
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
