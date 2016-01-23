// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
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

#include "dashmm.h"

#include "user_expansion.h"


struct UserSourceData {
  double pos[3];
  double mass;
};

struct UserTargetData {
  double pos[3];
  double phi[2];    //real, imag
};


struct InputArguments {
  int source_count;
  int target_count;
  int refinement_limit;
  int accuracy;
};


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
  fprintf(stdout, "Testing User expansion:\n");
  fprintf(stdout, "%d sources\n", retval.source_count);
  fprintf(stdout, "%d targets\n", retval.target_count);
  fprintf(stdout, "threshold: %d\n\n", retval.refinement_limit);

  return 0;
}


void pick_cube_position(double *pos) {
  pos[0] = (double)rand() / RAND_MAX;
  pos[1] = (double)rand() / RAND_MAX;
  pos[2] = (double)rand() / RAND_MAX;
}

double pick_mass(bool use_negative) {
  double retval = (double)rand() / RAND_MAX + 1.0;
  if (use_negative && (rand() % 2)) {
    retval *= -1.0;
  }
  return retval;
}


void set_sources(UserSourceData *sources, int source_count) {
  for (int i = 0; i < source_count; ++i) {
    pick_cube_position(sources[i].pos);
    sources[i].mass = pick_mass(true);
  }
}


void set_targets(UserTargetData *targets, int target_count) {
  for (int i = 0; i < target_count; ++i) {
    pick_cube_position(targets[i].pos);
    targets[i].phi[0] = 0.0;
    targets[i].phi[1] = 0.0;
  }
}


void perform_evaluation_test(InputArguments args) {
  srand(123456);

  //register User with DASHMM
  register_user_with_dashmm();

  //create some arrays
  UserSourceData *sources = static_cast<UserSourceData *>(
        malloc(sizeof(UserSourceData) * args.source_count));
  UserTargetData *targets = static_cast<UserTargetData *>(
        malloc(sizeof(UserTargetData) * args.target_count));

  set_sources(sources, args.source_count);
  set_targets(targets, args.target_count);

  //prep sources
  dashmm::ObjectHandle source_handle;
  auto err = dashmm::allocate_array(args.source_count, sizeof(UserSourceData),
                            &source_handle);
  assert(err == dashmm::kSuccess);
  err = dashmm::array_put(source_handle, 0, args.source_count, sources);
  assert(err == dashmm::kSuccess);

  //prep targets
  dashmm::ObjectHandle target_handle;
  err = dashmm::allocate_array(args.target_count, sizeof(UserTargetData),
                               &target_handle);
  assert(err == dashmm::kSuccess);
  err = dashmm::array_put(target_handle, 0, args.target_count, targets);
  assert(err == dashmm::kSuccess);

  //get method and expansion
  auto method = dashmm::fmm_method();
  dashmm::Expansion *test_expansion{
    new User{dashmm::Point{0.0, 0.0, 0.0}, args.accuracy}
  };

  assert(method && test_expansion);


  err = dashmm::evaluate(source_handle, offsetof(UserSourceData, pos),
                         offsetof(UserSourceData, mass),
                         target_handle, offsetof(UserTargetData, pos),
                         offsetof(UserTargetData, phi),
                         args.refinement_limit,
                         std::unique_ptr<dashmm::Method>{method},
                         std::unique_ptr<dashmm::Expansion>{test_expansion});

  //free up resources
  err = dashmm::deallocate_array(source_handle);
  assert(err == dashmm::kSuccess);
  err = dashmm::deallocate_array(target_handle);
  assert(err == dashmm::kSuccess);

  free(sources);
  free(targets);
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
