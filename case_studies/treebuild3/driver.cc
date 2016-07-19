#include <iostream>
#include <getopt.h>
#include "tree.h"

void usage(std::string exec) {
  std::cout << "Usage: " << exec << "--datatype=[c/s] " 
            << "--nsrc-per-rank=num " << "--ntar-per-rank=num "
            << "--threshold=num " << std::endl; 
}

int main(int argc, char *argv[]) {
  // default values 
  char datatype = 'c';
  int nsrc_per_rank = 1; 
  int ntar_per_rank = 1; 
  int threshold = 1; 
 
  if (hpx_init(&argc, &argv)) {
    std::cout << "HPX: failed to initialize" << std::endl;
  }

  // Parse command line 
  int opt = 0; 
  static struct option long_options[] = {
    {"datatype", required_argument, 0, 'd'}, 
    {"nsrc-per-rank", required_argument, 0, 's'}, 
    {"ntar-per-rank", required_argument, 0, 't'}, 
    {"threshold", required_argument, 0, 'l'}, 
    {"help", no_argument, 0, 'h'}, 
    {0, 0, 0, 0}
  }; 

  int long_index = 0; 
  while ((opt = getopt_long(argc, argv, "s:t:l:d:h", 
                            long_options, &long_index)) != -1) {
    switch (opt) {
    case 's':
      nsrc_per_rank = atoi(optarg); 
      break;
    case 't':
      ntar_per_rank = atoi(optarg); 
      break;
    case 'l':
      threshold = atoi(optarg); 
      break;
    case 'd':
      datatype = *optarg;
      break;
    case 'h':
      usage(std::string{argv[0]});
      return -1;
    }
  }

  if (hpx_run(&main_action, &nsrc_per_rank, &ntar_per_rank, 
              &datatype, &threshold)) {
    std::cout << "Failed to run main action" << std::endl;
  }

  hpx_finalize(); 
  return 0;
}
