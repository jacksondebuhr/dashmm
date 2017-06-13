//Code borrowed from stepping.cc and testmain.cc
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <getopt.h>
#include <sys/time.h>

#include <algorithm>
#include <map>
#include <memory>
#include <string>

#include "dashmm/dashmm.h"
struct Data{
	dashmm::Point position;
	double charge;
	std::complex<double> phi;
	int index;
};

struct InputArguments{
  int count;
  int refinement_limit;
  int steps;
  std::string output;
};

//Create evaluator to use the same source and target points
dashmm::Evaluator<Data, Data, dashmm::LaplaceCOM, dashmm::BH> bheval{};

//Parse command line inputs (stepping.cc)
int read_arguments(int argc, char **argv, InputArguments &retval) {
  // Set defaults
  retval.count = 10000;
  retval.refinement_limit = 40;
  retval.steps = 20;
  retval.output.clear();

  int opt = 0;
  static struct option long_options[] = {
    {"nsources", required_argument, 0, 's'},
    {"threshold", required_argument, 0, 'l'},
    {"nsteps", required_argument, 0, 'p'},
    {"output", required_argument, 0, 'o'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
  };

  int long_index = 0;
  while ((opt = getopt_long(argc, argv, "m:s:w:t:g:l:v:a:h",
                            long_options, &long_index)) != -1) {
    std::string verifyarg{ };
    switch (opt) {
    case 's':
      retval.count = atoi(optarg);
      break;
    case 'l':
      retval.refinement_limit = atoi(optarg);
      break;
    case 'p':
      retval.steps = atoi(optarg);
      break;
    case 'o':
      retval.output = std::string(optarg);
      break;
    case 'h':
      print_usage(argv[0]);
      return -1;
    case '?':
      return -1;
    }
  }

  // test the inputs
  if (retval.count < 1) {
    fprintf(stderr, "Usage ERROR: nsources must be positive.\n");
    return -1;
  }

  // print out summary
  if (dashmm::get_my_rank() == 0) {
    fprintf(stdout, "Testing DASHMM:\n");
    fprintf(stdout, "%d sources taking %d steps\n", retval.count, retval.steps);
    fprintf(stdout, "threshold: %d\n", retval.refinement_limit);
    if (!retval.output.empty()) {
      fprintf(stdout, "output in file: %s\n\n", retval.output.c_str());
    }
  } else {
    retval.count = 0;
  }

  return 0;
}


//The main driver: 
void perform_evaluation_test(InputArguments args){

}


int main(int argc, char** argv){
	auto err = dashmm::init(&argc, &argv);
	assert(err == dashmm::kSuccess);

	InputArguments inputargs;
	int usage_error = read_arguments(argc, argv, inputargs);

	if(!usage_error){
		perform_evaluation_test(inputargs);
	}

	err = dashmm::finalize();
	assert(err == dashmm::kSuccess);

	return 0;
}
