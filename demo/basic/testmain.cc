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
#include <complex>
#include <map>
#include <memory>
#include <string>

#include "dashmm/dashmm.h"


struct SourceData {
  dashmm::Point position;
  double charge;
};

struct TargetData {
  dashmm::Point position;
  std::complex<double> phi;    //real, imag
  int index;
};


// Here we create the three evaluator objects that we shall need in this demo.
// These must be instantiated before the call to dashmm::init so that they
// might register the relevant actions with the runtime system.
dashmm::Evaluator<SourceData, TargetData,
                  dashmm::LaplaceCOM, dashmm::BH> bheval{};
dashmm::Evaluator<SourceData, TargetData,
                  dashmm::Laplace, dashmm::FMM> fmmeval{};
dashmm::Evaluator<SourceData, TargetData, 
                  dashmm::Laplace, dashmm::FMM97> fmm97eval{}; 
dashmm::Evaluator<SourceData, TargetData,
                  dashmm::LaplaceCOM, dashmm::Direct> directeval{};



struct InputArguments {
  int source_count;
  std::string source_type;
  int target_count;
  std::string target_type;
  int refinement_limit;
  std::string test_case;
  bool verify;
  int accuracy;
};


void print_usage(char *progname) {
  fprintf(stdout, "Usage: %s [OPTIONS]\n\n"
          "Options available: [possible/values] (default value)\n"
          "  --method=[fmm/fmm97/bh]            method to use (fmm)\n"
          "  --nsources=num               "
          "number of source points to generate (10000)\n"
          "  --sourcedata=[cube/sphere/plummer]\n"
          "                               source distribution type (cube)\n"
          "  --ntargets=num               "
          "number of target points to generate (10000)\n"
          "  --targetdata=[cube/sphere/plummer]\n"
          "                               target distribution type (cube)\n"
          "  --threshold=num              "
          "source and target tree partition refinement\n"
          "                                 limit (40)\n"
          "  --accuracy=num               "
          "number of digits of accuracy for fmm (3)\n"
          "  --verify=[yes/no]            "
          "perform an accuracy test comparing to direct\n"
          "                                 summation (yes)\n", progname);
}

int read_arguments(int argc, char **argv, InputArguments &retval) {
  //Set defaults
  retval.source_count = 10000;
  retval.source_type = std::string{"cube"};
  retval.target_count = 10000;
  retval.target_type = std::string{"cube"};
  retval.refinement_limit = 40;
  retval.test_case = std::string{"fmm97"};
  retval.verify = true;
  retval.accuracy = 3;

  int opt = 0;
  static struct option long_options[] = {
    {"method", required_argument, 0, 'm'},
    {"nsources", required_argument, 0, 's'},
    {"sourcedata", required_argument, 0, 'w'},
    {"ntargets", required_argument, 0, 't'},
    {"targetdata", required_argument, 0, 'g'},
    {"threshold", required_argument, 0, 'l'},
    {"verify", required_argument, 0, 'v'},
    {"accuracy", required_argument, 0, 'a'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
  };

  int long_index = 0;
  while ((opt = getopt_long(argc, argv, "m:s:w:t:g:l:v:a:h",
                            long_options, &long_index)) != -1) {
    std::string verifyarg{};
    switch (opt) {
    case 'm':
      retval.test_case = optarg;
      break;
    case 's':
      retval.source_count = atoi(optarg);
      break;
    case 'w':
      retval.source_type = optarg;
      break;
    case 't':
      retval.target_count = atoi(optarg);
      break;
    case 'g':
      retval.target_type = optarg;
      break;
    case 'l':
      retval.refinement_limit = atoi(optarg);
      break;
    case 'v':
      verifyarg = optarg;
      retval.verify = (verifyarg == std::string{"yes"});
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
  if (retval.test_case != "bh" && retval.test_case != "fmm" 
      && retval.test_case != "fmm97") {
    fprintf(stderr, "Usage ERROR: unknown method '%s'\n",
            retval.test_case.c_str());
    return -1;
  }
  if (retval.source_type != "cube" && retval.source_type != "sphere"
        && retval.source_type != "plummer") {
    fprintf(stderr, "Usage ERROR: unknown source type '%s'\n",
            retval.source_type.c_str());
    return -1;
  }
  if (retval.target_type != "cube" && retval.target_type != "sphere"
        && retval.target_type != "plummer") {
    fprintf(stderr, "Usage ERROR: unknown target type '%s'\n",
            retval.target_type.c_str());
    return -1;
  }

  //print out summary
  if (hpx_get_my_rank() == 0) {
    fprintf(stdout, "Testing DASHMM:\n");
    fprintf(stdout, "%d sources in a %s distribution\n",
            retval.source_count, retval.source_type.c_str());
    fprintf(stdout, "%d targets in a %s distribution\n",
            retval.target_count, retval.target_type.c_str());
    fprintf(stdout, "method: %s \nthreshold: %d\n\n",
            retval.test_case.c_str(), retval.refinement_limit);
  } else {
    // Only have rank 0 create data
    retval.source_count = 0;
    retval.target_count = 0;
  }

  return 0;
}


inline double getticks(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double) (tv.tv_sec * 1e6 + tv.tv_usec);
}

inline double elapsed(double t1, double t0) {
  return (double) (t1 - t0);
}


dashmm::Point pick_cube_position() {
  double pos[3];
  pos[0] = (double)rand() / RAND_MAX - 0.5;
  pos[1] = (double)rand() / RAND_MAX - 0.5;
  pos[2] = (double)rand() / RAND_MAX - 0.5;
  return dashmm::Point{pos[0], pos[1], pos[2]};
}


dashmm::Point pick_sphere_position() {
  double r = 1.0;
  double ctheta = 2.0 * (double)rand() / RAND_MAX - 1.0;
  double stheta = sqrt(1.0 - ctheta * ctheta);
  double phi = 2.0 * 3.1415926535 * (double)rand() / RAND_MAX;
  double pos[3];
  pos[0] = r * stheta * cos(phi);
  pos[1] = r * stheta * sin(phi);
  pos[2] = r * ctheta;
  return dashmm::Point{pos[0], pos[1], pos[2]};
}


double pick_mass(bool use_negative) {
  double retval = (double)rand() / RAND_MAX + 1.0;
  if (use_negative && (rand() % 2)) {
    retval *= -1.0;
  }
  //double retval = 1.0 * rand() / RAND_MAX - 0.5; 
  return retval;
}


dashmm::Point pick_plummer_position() {
  //NOTE: This is using a = 1
  double unif = (double)rand() / RAND_MAX;
  double r = 1.0 / sqrt(pow(unif, -2.0 / 3.0) - 1);
  double ctheta = 2.0 * (double)rand() / RAND_MAX - 1.0;
  double stheta = sqrt(1.0 - ctheta * ctheta);
  double phi = 2.0 * 3.1415926535 * (double)rand() / RAND_MAX;
  double pos[3];
  pos[0] = r * stheta * cos(phi);
  pos[1] = r * stheta * sin(phi);
  pos[2] = r * ctheta;
  return dashmm::Point{pos[0], pos[1], pos[2]};
}


double pick_plummer_mass(int count) {
  //We take the total mass to be 100
  return 100.0 / count;
}


void set_sources(SourceData *sources, int source_count,
                 std::string source_type, std::string test_case) {
  bool use_negative = (test_case != std::string{"bh"}); 
  if (source_type == std::string{"cube"}) {
    for (int i = 0; i < source_count; ++i) {
      sources[i].position = pick_cube_position();
      sources[i].charge = pick_mass(use_negative);
    }
  } else if (source_type == std::string{"sphere"}) {
    //Sphere
    for (int i = 0; i < source_count; ++i) {
      sources[i].position = pick_sphere_position();
      sources[i].charge = pick_mass(use_negative);
    }
  } else {
    //Plummer
    for (int i = 0; i < source_count; ++i) {
      sources[i].position = pick_plummer_position();
      sources[i].charge = pick_plummer_mass(source_count);
    }
  }
}


void set_targets(TargetData *targets, int target_count,
                 std::string target_type) {
  if (target_type == std::string{"cube"}) {
    //Cube
    for (int i = 0; i < target_count; ++i) {
      targets[i].position = pick_cube_position();
      targets[i].phi = 0.0;
      targets[i].index = i;
    }
  } else if (target_type == std::string{"sphere"}) {
    //Sphere
    for (int i = 0; i < target_count; ++i) {
      targets[i].position = pick_sphere_position();
      targets[i].phi = 0.0;
      targets[i].index = i;
    }
  } else {
    //Plummer
    for (int i = 0; i < target_count; ++i) {
      targets[i].position = pick_plummer_position();
      targets[i].phi = 0.0;
      targets[i].index = i;
    }
  }
}


void compare_results(TargetData *targets, int target_count,
                     TargetData *exacts, int exact_count) {
  if (hpx_get_my_rank()) return;

  //create a map from index into offset fort targets
  std::map<int, int> offsets{};
  for (int i = 0; i < target_count; ++i) {
    offsets[targets[i].index] = i;
  }

  //now we loop over the exact results and compare
  double numerator = 0.0;
  double denominator = 0.0;
  double maxrel = 0.0;
  for (int i = 0; i < exact_count; ++i) {
    auto j = offsets.find(exacts[i].index);
    assert(j != offsets.end());
    int idx = j->second;
    double relerr = fabs(targets[idx].phi.real() - exacts[i].phi.real());
    numerator += relerr * relerr;
    denominator += exacts[i].phi.real() * exacts[i].phi.real();
    if (relerr / exacts[i].phi.real() > maxrel) {
      maxrel = relerr / exacts[i].phi.real();
    }
  }
  fprintf(stdout, "Error for %d test points: %4.3e (max %4.3e)\n",
                  exact_count, sqrt(numerator / denominator), maxrel);
}


void perform_evaluation_test(InputArguments args) {
  srand(123456);

  //create some arrays
  SourceData *sources{nullptr};
  TargetData *targets{nullptr};
  if (args.source_count) {
    sources = reinterpret_cast<SourceData *>(
        new char [sizeof(SourceData) * args.source_count]);
    set_sources(sources, args.source_count, args.source_type, args.test_case);
  }
  if (args.target_count) {
    targets = reinterpret_cast<TargetData *>(
        new char [sizeof(TargetData) * args.target_count]);
    set_targets(targets, args.target_count, args.target_type);
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

  //save a few targets in case of direct comparison
  int test_count{0};
  if (hpx_get_my_rank() == 0) {
    test_count = 400;
  }
  if (test_count > args.target_count) {
    test_count = args.target_count;
  }
  TargetData *test_targets{nullptr};
  if (test_count) {
    test_targets = reinterpret_cast<TargetData *>(
          new char [sizeof(TargetData) * test_count]);

    for (int i = 0; i < test_count; ++i) {
      int idx = i * (args.target_count / test_count);
      assert(idx < args.target_count);
      test_targets[i] = targets[idx];
    }
  }

  //Perform the evaluation
  double t0{};
  double tf{};
  if (args.test_case == std::string{"bh"}) {
    dashmm::LaplaceCOM<SourceData, TargetData> expansion{
        dashmm::Point{0.0, 0.0, 0.0}, 0, dashmm::kNoRoleNeeded};
    dashmm::BH<SourceData, TargetData, dashmm::LaplaceCOM> method{0.6};

    t0 = getticks();
    err = bheval.evaluate(source_handle, target_handle, args.refinement_limit,
                          method, expansion);
    assert(err == dashmm::kSuccess);
    tf = getticks();
  } else if (args.test_case == std::string{"fmm"}) {
    dashmm::Laplace<SourceData, TargetData> expansion{
          dashmm::Point{0.0, 0.0, 0.0}, args.accuracy, dashmm::kNoRoleNeeded};
    dashmm::FMM<SourceData, TargetData, dashmm::Laplace> method{};

    t0 = getticks();
    err = fmmeval.evaluate(source_handle, target_handle, args.refinement_limit,
                           method, expansion);
    assert(err == dashmm::kSuccess);
    tf = getticks();
  } else if (args.test_case == std::string{"fmm97"}) {
    dashmm::Laplace<SourceData, TargetData> expansion{
      dashmm::Point{0.0, 0.0, 0.0}, args.accuracy, dashmm::kNoRoleNeeded}; 
    dashmm::FMM97<SourceData, TargetData, dashmm::Laplace> method{}; 

    t0 = getticks(); 
    err = fmm97eval.evaluate(source_handle, target_handle, 
                             args.refinement_limit, method, expansion); 
    assert(err == dashmm::kSuccess); 
    tf = getticks(); 
  }

  fprintf(stdout, "Evaluation took %lg [us]\n", elapsed(tf, t0));

  if (args.verify) {
    // Create array for test targets
    dashmm::Array<TargetData> test_handle{};
    err = test_handle.allocate(test_count);
    assert(err == dashmm::kSuccess);
    err = test_handle.put(0, test_count, test_targets);
    assert(err == dashmm::kSuccess);

    //do direct evaluation
    dashmm::Direct<SourceData, TargetData, dashmm::LaplaceCOM> direct{};
    dashmm::LaplaceCOM<SourceData, TargetData> direxp{
        dashmm::Point{0.0, 0.0, 0.0}, 0, dashmm::kNoRoleNeeded};
    err = directeval.evaluate(source_handle, test_handle, args.refinement_limit,
                              direct, direxp);
    assert(err == dashmm::kSuccess);

    //Get the results from the global address space
    err = target_handle.get(0, args.target_count, targets);
    assert(err == dashmm::kSuccess);
    err = test_handle.get(0, test_count, test_targets);
    assert(err == dashmm::kSuccess);

    //Test error
    compare_results(targets, args.target_count, test_targets, test_count);

    err = test_handle.destroy();
    assert(err == dashmm::kSuccess);
    delete [] test_targets;
  }

  //free up resources
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
