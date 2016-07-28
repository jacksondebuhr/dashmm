#include <iostream>
#include <getopt.h>
#include "tree.h"

void usage(std::string exec) {
  std::cout << "Usage: " << exec << " --scaling=[w/s] " 
            << "--datatype=[c/s] " << " --exchange=[y/n] "
            << "--nsrc=num " << "--ntar=num " 
            << "--threshold=num " << "--nseed=num " << std::endl; 
  std::cout << "Note: --nseed=num required if --scaling=s" << std::endl;
}

int main(int argc, char *argv[]) {
  // default values 
  char scaling = 'w'; // strong scale, 'w' for weak scaling
  char datatype = 'c'; // cubic point distribution, 's' for sphere distrbution
  char exchange = 'n'; // exchange trees
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
    {"exchange", required_argument, 0, 'e'}, 
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
    case 'e':
      exchange = *optarg;
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
    if (hpx_run(&main_action, &scaling, &datatype, &exchange, 
                &nsrc, &ntar, &threshold, &nseed)) {
      std::cout << "Failed to run main action" << std::endl;
    }
  }

  hpx_finalize(); 
  return 0;
}
