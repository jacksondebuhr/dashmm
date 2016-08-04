#include <iostream>
#include <getopt.h>
#include "tree.h"

void usage(std::string exec) {
  std::cout << "Usage: " << exec << " --scaling=[w/s] "
            << "--datatype=[c/s] "
            << "--nsrc=num " << "--ntar=num "
            << "--threshold=num " << "--nseed=num " << std::endl;
  std::cout << "Note: --nseed=num required if --scaling=s" << std::endl;
}

/// This will generate the points for the trees.
/// This is pretty simple. The local offset to the correct part of the
/// meta data for the array is computed, and the values are filled in.
/// nothing fancy.
int generate_input_handler(char scaling, char datatype, int nsrc, int ntar,
                           int nseed, hpx_addr_t sources, hpx_addr_t targets) {
  int rank = hpx_get_my_rank();
  int num_ranks = hpx_get_num_ranks();
  hpx_addr_t curr_s = hpx_addr_add(sources, sizeof(ArrayMetaData) * rank,
                                   sizeof(ArrayMetaData));
  hpx_addr_t curr_t = hpx_addr_add(targets, sizeof(ArrayMetaData) * rank,
                                   sizeof(ArrayMetaData));

  ArrayMetaData *meta_s{nullptr}, *meta_t{nullptr};
  if (!hpx_gas_try_pin(curr_s, (void **)&meta_s) ||
      !hpx_gas_try_pin(curr_t, (void **)&meta_t))
    return HPX_ERROR;

  if (scaling == 'w') {
    // Generate input for weak scaling test
    // nsrc and ntar are the number of sources and targets per rank

    // Set the seed of the random number generator. Note that if the seed is set
    // to 1, the generator is reinitialized to its initial value. For this
    // reason, the seed is chosen to be 2 + rank here.
    generate_weak_scaling_input(meta_s, nsrc, meta_t, ntar, datatype, rank + 2);
  } else {
    // Generate input for strong scaling test
    // nsrc and ntar are the total number of sources and targets
    generate_strong_scaling_input(meta_s, nsrc, meta_t, ntar, datatype,
                                  rank, num_ranks, nseed);
  }
  hpx_gas_unpin(curr_s);
  hpx_gas_unpin(curr_t);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, generate_input_action, generate_input_handler,
           HPX_CHAR, HPX_CHAR, HPX_INT, HPX_INT, HPX_INT, HPX_ADDR, HPX_ADDR);

/// Pretty simple really, just gets rid of the sources and targets
int destroy_input_handler(hpx_addr_t sources, hpx_addr_t targets) {
  int rank = hpx_get_my_rank();
  hpx_addr_t curr_s = hpx_addr_add(sources, sizeof(ArrayMetaData) * rank,
                                   sizeof(ArrayMetaData));
  hpx_addr_t curr_t = hpx_addr_add(targets, sizeof(ArrayMetaData) * rank,
                                   sizeof(ArrayMetaData));

  ArrayMetaData *meta_s{nullptr}, *meta_t{nullptr};
  if (!hpx_gas_try_pin(curr_s, (void **)&meta_s) ||
      !hpx_gas_try_pin(curr_t, (void **)&meta_t))
    return HPX_ERROR;

  delete [] meta_s->data;
  delete [] meta_t->data;

  hpx_gas_unpin(curr_s);
  hpx_gas_unpin(curr_t);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, destroy_input_action, destroy_input_handler,
           HPX_ADDR, HPX_ADDR);

/// The main action of the code starts here
int main_handler(char scaling, char datatype,
                 int nsrc, int ntar, int threshold, int nseed) {
  int num_ranks = hpx_get_num_ranks();

  // Generate input data
  hpx_addr_t sources = hpx_gas_calloc_cyclic(num_ranks,
                                             sizeof(ArrayMetaData), 0);
  hpx_addr_t targets = hpx_gas_calloc_cyclic(num_ranks,
                                             sizeof(ArrayMetaData), 0);

  hpx_bcast_rsync(generate_input_action, &scaling, &datatype, &nsrc,
                  &ntar, &nseed, &sources, &targets);

  // Partition points to create dual tree

  hpx_time_t timer_start = hpx_time_now();

  // Determine domain geometry
  hpx_addr_t domain_geometry =
    hpx_lco_reduce_new(num_ranks, sizeof(double) * 6,
                       domain_geometry_init_action,
                       domain_geometry_op_action);

  hpx_bcast_lsync(set_domain_geometry_action, HPX_NULL,
                  &sources, &targets, &domain_geometry);

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
  unif_done = hpx_gas_calloc_cyclic(num_ranks, sizeof(hpx_addr_t), 0);
  unif_grid = hpx_gas_calloc_cyclic(num_ranks, sizeof(ArrayMetaData), 0);
  sorted_src = hpx_gas_calloc_cyclic(num_ranks, sizeof(ArrayMetaData), 0);
  sorted_tar = hpx_gas_calloc_cyclic(num_ranks, sizeof(ArrayMetaData), 0);

  hpx_bcast_rsync(init_partition_action, &unif_count, &unif_done, &unif_grid,
                  &sorted_src, &sorted_tar, &unif_level, &threshold,
                  &corner_x, &corner_y, &corner_z, &size);

  // Measure the partitioning setup. This is all just getting this and that
  // ready to do. No real work is done in the above.
  hpx_time_t timer_middle = hpx_time_now();

  hpx_bcast_rsync(create_dual_tree_action, &sources, &targets);

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

  // Destroy input data
  hpx_bcast_rsync(destroy_input_action, &sources, &targets);

  hpx_lco_delete_sync(unif_count);
  hpx_gas_free_sync(unif_done);
  hpx_gas_free_sync(unif_grid);
  hpx_gas_free_sync(sorted_src);
  hpx_gas_free_sync(sorted_tar);
  hpx_gas_free_sync(sources);
  hpx_gas_free_sync(targets);

  hpx_exit(0, nullptr);
}
HPX_ACTION(HPX_DEFAULT, 0, main_action, main_handler,
           HPX_CHAR, HPX_CHAR, HPX_INT, HPX_INT, HPX_INT, HPX_INT);

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

  if (valid_arguments) {
    if (hpx_run(&main_action, nullptr, &scaling, &datatype,
                &nsrc, &ntar, &threshold, &nseed)) {
      std::cout << "Failed to run main action" << std::endl;
    }
  }

  hpx_finalize();
  return 0;
}
