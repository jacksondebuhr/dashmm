// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2017, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license. See the LICENSE file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================


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

#include "user_expansion.h"


// Officially speaking, the User expansion type does not place any requirements
// on Source types or Target types (see user_expansion.h). However, we generate
// positions and masses for the data anyway.

struct SourceData {
  dashmm::Point position;
  double mass;
};

struct TargetData {
  dashmm::Point position;    //real, imag
  dashmm::dcomplex_t phi;
};


// These are the input arguments to the program that are configurable on the
// command line.
struct InputArguments {
  int source_count;
  int target_count;
  int refinement_limit;
  int accuracy;
};


// To use an instance of DASHMM for a set of types, one must create one (and
// only one) instance of the assicated Evaluator type. This object must be
// created prior to dashmm::init() being called.
dashmm::Evaluator<SourceData, TargetData, User, dashmm::FMM> usereval{};


void print_usage(char *progname) {
  fprintf(stdout, "Usage: %s [OPTIONS]\n\n"
"Options available: [possible/values] (default value)\n"
"  --nsources=num               number of source points to generate (10)\n"
"  --ntargets=num               number of target points to generate (10)\n"
"  --threshold=num              source and target tree partition refinement\n"
"                                 limit (1)\n"
"  --accuracy=num               number of digits of accuracy (3)\n",
    progname);
}


int read_arguments(int argc, char **argv, InputArguments &retval) {
  //Set defaults
  retval.source_count = 10;
  retval.target_count = 10;
  retval.refinement_limit = 1;
  retval.accuracy = 3;

  int opt = 0;
  static struct option long_options[] = {
    {"nsources", required_argument, 0, 's'},
    {"ntargets", required_argument, 0, 't'},
    {"threshold", required_argument, 0, 'l'},
    {"accuracy", required_argument, 0, 'a'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
  };

  int long_index = 0;
  while ((opt = getopt_long(argc, argv, "m:s:w:t:g:l:v:a:h",
                            long_options, &long_index)) != -1) {
    std::string verifyarg{};
    switch (opt) {
    case 's':
      retval.source_count = atoi(optarg);
      break;
    case 't':
      retval.target_count = atoi(optarg);
      break;
    case 'l':
      retval.refinement_limit = atoi(optarg);
      break;
    case 'a':
      retval.accuracy = atoi(optarg);
      break;
    case 'h':
      print_usage(argv[0]);
      return -1;
    case '?':
      return -1;
    }
  }

  //test the inputs
  if (retval.source_count < 1) {
    fprintf(stderr, "Usage ERROR: nsources must be positive.\n");
    return -1;
  }
  if (retval.target_count < 1) {
    fprintf(stderr, "Usage ERROR: ntargets must be positive\n");
    return -1;
  }

  //print out summary
  if (dashmm::get_my_rank() == 0) {
    fprintf(stdout, "Testing User expansion:\n");
    fprintf(stdout, "%d sources\n", retval.source_count);
    fprintf(stdout, "%d targets\n", retval.target_count);
    fprintf(stdout, "threshold: %d\n\n", retval.refinement_limit);
  } else {
    retval.source_count = 0;
    retval.target_count = 0;
  }

  return 0;
}


dashmm::Point pick_cube_position() {
  double pos[3];
  pos[0] = (double)rand() / RAND_MAX;
  pos[1] = (double)rand() / RAND_MAX;
  pos[2] = (double)rand() / RAND_MAX;
  return dashmm::Point{pos[0], pos[1], pos[2]};
}

double pick_mass(bool use_negative) {
  double retval = (double)rand() / RAND_MAX + 1.0;
  if (use_negative && (rand() % 2)) {
    retval *= -1.0;
  }
  return retval;
}


void set_sources(SourceData *sources, int source_count) {
  for (int i = 0; i < source_count; ++i) {
    sources[i].position = pick_cube_position();
    sources[i].mass = pick_mass(true);
  }
}


void set_targets(TargetData *targets, int target_count) {
  for (int i = 0; i < target_count; ++i) {
    targets[i].position = pick_cube_position();
    targets[i].phi = 0.0;
  }
}


void perform_evaluation_test(InputArguments args) {
  srand(123456);

  //create some arrays
  SourceData *sources{nullptr};
  if (args.source_count) {
    sources = reinterpret_cast<SourceData *>(
          new char [sizeof(SourceData) * args.source_count]);
    set_sources(sources, args.source_count);
  }
  TargetData *targets{nullptr};
  if (args.target_count) {
    targets = reinterpret_cast<TargetData *>(
          new char [sizeof(TargetData) * args.target_count]);
    set_targets(targets, args.target_count);
  }

  //prep sources
  dashmm::Array<SourceData> source_handle{};
  int err = source_handle.allocate(args.source_count);
  assert(err == dashmm::kSuccess);
  err = source_handle.put(0, args.source_count, sources);
  assert(err == dashmm::kSuccess);

  //prep targets
  dashmm::Array<TargetData> target_handle{};
  err = target_handle.allocate(args.target_count);
  assert(err == dashmm::kSuccess);
  err = target_handle.put(0, args.target_count, targets);
  assert(err == dashmm::kSuccess);

  // Create method - for this example we use FMM to showcase more operations.
  // Methods are explicitly passed into evaluate() so that any parameters of
  // the method might propagate through the evaluation.
  dashmm::FMM<SourceData, TargetData, User> method{};

  // All that is left is to call evaluate from the Evaluator object.
  std::vector<double> kparm{};
  err = usereval.evaluate(source_handle, target_handle, args.refinement_limit,
                          &method, args.accuracy, &kparm);

  // Clean up resources
  err = source_handle.destroy();
  assert(err == dashmm::kSuccess);
  err = target_handle.destroy();
  assert(err == dashmm::kSuccess);

  delete [] sources;
  delete [] targets;
}


int main(int argc, char **argv) {
  auto err = dashmm::init(&argc, &argv);
  assert(err == dashmm::kSuccess);

  InputArguments inputargs;
  int usage_error = read_arguments(argc, argv, inputargs);

  if (!usage_error) {
    perform_evaluation_test(inputargs);
  }

  err = dashmm::finalize();
  assert(err == dashmm::kSuccess);

  return 0;
}
