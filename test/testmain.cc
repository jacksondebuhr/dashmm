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
  std::string source_type;
  int target_count;
  std::string target_type;
  int refinement_limit;
  std::string test_case;
  bool verify;
  int accuracy;
};


void print_usage(char *progname) {
  //TODO: Improve formatting; add other methods when available
  fprintf(stdout, "Usage: %s [OPTIONS]\n\n"
"Options available: [possible/values] (default value)\n"
"  --method=[fmm/bh]            method to use (fmm)\n"
"  --nsources=num               number of source points to generate (10000)\n"
"  --sourcedata=[cube/sphere/plummer]\n"
"                               source distribution type (cube)\n"
"  --ntargets=num               number of target points to generate (10000)\n"
"  --targetdata=[cube/sphere/plummer]\n"
"                               target distribution type (cube)\n"
"  --threshold=num              source and target tree partition refinement\n"
"                                 limit (40)\n"
"  --accuracy=num               number of digits of accuracy for fmm (3)\n"
"  --verify=[yes/no]            perform an accuracy test comparing to direct\n"
"                                 summation (yes)\n", progname);
}


int read_arguments(int argc, char **argv, InputArguments &retval) {
  //Set defaults
  retval.source_count = 10000;
  retval.source_type = std::string{"cube"};
  retval.target_count = 10000;
  retval.target_type = std::string{"cube"};
  retval.refinement_limit = 40;
  retval.test_case = std::string{"fmm"};
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
  //TODO add other cases when available
  if (retval.test_case != "bh" && retval.test_case != "fmm") {
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
  fprintf(stdout, "Testing DASHMM:\n");
  fprintf(stdout, "%d sources in a %s distribution\n",
          retval.source_count, retval.source_type.c_str());
  fprintf(stdout, "%d targets in a %s distribution\n",
          retval.target_count, retval.target_type.c_str());
  fprintf(stdout, "method: %s \nthreshold: %d\n\n",
          retval.test_case.c_str(), retval.refinement_limit);

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


void pick_cube_position(double *pos) {
  pos[0] = (double)rand() / RAND_MAX;
  pos[1] = (double)rand() / RAND_MAX;
  pos[2] = (double)rand() / RAND_MAX;
}


void pick_sphere_position(double *pos) {
  double r = 1.0;
  double ctheta = 2.0 * (double)rand() / RAND_MAX - 1.0;
  double stheta = sqrt(1.0 - ctheta * ctheta);
  double phi = 2.0 * 3.1415926535 * (double)rand() / RAND_MAX;
  pos[0] = r * stheta * cos(phi);
  pos[1] = r * stheta * sin(phi);
  pos[2] = r * ctheta;
}


double pick_mass() {
  return (double)rand() / RAND_MAX + 1.0;
}


void pick_plummer_position(double *pos) {
  //NOTE: This is using a = 1
  double unif = (double)rand() / RAND_MAX;
  double r = 1.0 / sqrt(pow(unif, -2.0 / 3.0) - 1);
  double ctheta = 2.0 * (double)rand() / RAND_MAX - 1.0;
  double stheta = sqrt(1.0 - ctheta * ctheta);
  double phi = 2.0 * 3.1415926535 * (double)rand() / RAND_MAX;
  pos[0] = r * stheta * cos(phi);
  pos[1] = r * stheta * sin(phi);
  pos[2] = r * ctheta;
}


double pick_plummer_mass(int count) {
  //We take the total mass to be 100
  return 100.0 / count;
}


void set_sources(UserSourceData *sources, int source_count,
                 std::string source_type) {
  if (source_type == std::string{"cube"}) {
    for (int i = 0; i < source_count; ++i) {
      pick_cube_position(sources[i].pos);
      sources[i].mass = pick_mass();
    }
  } else if (source_type == std::string{"sphere"}) {
    //Sphere
    for (int i = 0; i < source_count; ++i) {
      pick_sphere_position(sources[i].pos);
      sources[i].mass = pick_mass();
    }
  } else {
    //Plummer
    for (int i = 0; i < source_count; ++i) {
      pick_plummer_position(sources[i].pos);
      sources[i].mass = pick_plummer_mass(source_count);
    }
  }
}


void set_targets(UserTargetData *targets, int target_count,
                 std::string target_type) {
  if (target_type == std::string{"cube"}) {
    //Cube
    for (int i = 0; i < target_count; ++i) {
      pick_cube_position(targets[i].pos);
      targets[i].phi[0] = 0.0;
      targets[i].phi[1] = 0.0;
    }
  } else if (target_type == std::string{"sphere"}) {
    //Sphere
    for (int i = 0; i < target_count; ++i) {
      pick_sphere_position(targets[i].pos);
      targets[i].phi[0] = 0.0;
      targets[i].phi[1] = 0.0;
    }
  } else {
    //Plummer
    for (int i = 0; i < target_count; ++i) {
      pick_plummer_position(targets[i].pos);
      targets[i].phi[0] = 0.0;
      targets[i].phi[1] = 0.0;
    }
  }
}


void perform_evaluation_test(InputArguments args) {
  srand(123456);

  //create some arrays
  UserSourceData *sources = static_cast<UserSourceData *>(
        malloc(sizeof(UserSourceData) * args.source_count));
  UserTargetData *targets = static_cast<UserTargetData *>(
        malloc(sizeof(UserTargetData) * args.target_count));

  set_sources(sources, args.source_count, args.source_type);
  set_targets(targets, args.target_count, args.target_type);

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
  dashmm::Method *test_method{nullptr};
  dashmm::Expansion *test_expansion{nullptr};
  if (args.test_case == std::string{"bh"}) {
    test_method = dashmm::bh_method(0.6);
    test_expansion = dashmm::laplace_COM_expansion();
  } else if (args.test_case == std::string{"fmm"}) {
    test_method = dashmm::fmm_method();
    test_expansion = dashmm::laplace_sph_expansion(args.accuracy);
  }

  assert(test_method && test_expansion);

  //evaluate - first for the approximate version
  double t0 = getticks();
  err = dashmm::evaluate(source_handle, offsetof(UserSourceData, pos),
                         offsetof(UserSourceData, mass),
                         target_handle, offsetof(UserTargetData, pos),
                         offsetof(UserTargetData, phi),
                         args.refinement_limit,
                         std::unique_ptr<dashmm::Method>{test_method},
                         std::unique_ptr<dashmm::Expansion>{test_expansion});
  double tf = getticks();
  fprintf(stdout, "Evaluation took %lg [us]\n", elapsed(tf, t0));

  if (args.verify) {
    //work out the test count
    int test_count = 400;
    if (test_count > args.target_count) {
      test_count = args.target_count;
    }
    UserTargetData *test_targets = static_cast<UserTargetData *>(
          malloc(sizeof(UserTargetData) * test_count));

    //Save a few targets for direct comparison
    for (int i = 0; i < test_count; ++i) {
      int idx = i * (args.target_count / test_count);
      assert(idx < args.target_count);
      test_targets[i] = targets[idx];
    }
    dashmm::ObjectHandle test_handle;
    err = dashmm::allocate_array(test_count, sizeof(UserTargetData),
                                 &test_handle);
    assert(err == dashmm::kSuccess);
    err = dashmm::array_put(test_handle, 0, test_count, test_targets);
    assert(err == dashmm::kSuccess);

    //do direct evaluation
    auto direct = dashmm::direct_method();
    auto direxp = dashmm::laplace_COM_expansion();
    err = dashmm::evaluate(source_handle, offsetof(UserSourceData, pos),
                           offsetof(UserSourceData, mass),
                           test_handle, offsetof(UserTargetData, pos),
                           offsetof(UserTargetData, phi),
                           args.refinement_limit,
                           std::unique_ptr<dashmm::Method>{direct},
                           std::unique_ptr<dashmm::Expansion>{direxp});

    err = dashmm::array_get(target_handle, 0, args.target_count, targets);
    assert(err == dashmm::kSuccess);
    err = dashmm::array_get(test_handle, 0, test_count, test_targets);
    assert(err == dashmm::kSuccess);
    err = dashmm::deallocate_array(test_handle);
    assert(err == dashmm::kSuccess);

    //Test error
    double numerator = 0.0;
    double denominator = 0.0;
    double maxrel = 0.0;
    for (int i = 0; i < test_count; ++i) {
      int idx = i * (args.target_count / test_count);
      assert(idx < args.target_count);
      double relerr = fabs(targets[idx].phi[0] - test_targets[i].phi[0]);
      numerator += relerr * relerr;
      denominator += test_targets[i].phi[0] * test_targets[i].phi[0];
      if (relerr / test_targets[i].phi[0] > maxrel) {
        maxrel = relerr / test_targets[i].phi[0];
      }
    }
    fprintf(stdout, "Error for %d test points: %lg (max %lg)\n",
                    test_count, sqrt(numerator / denominator), maxrel);

    free(test_targets);
  }

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
