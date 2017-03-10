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
#include <complex>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include "dashmm/dashmm.h"


constexpr int kNSources = 10000;
constexpr int kNTargets = 10000;
constexpr int kRefinementLimit = 40;
constexpr double kYukawaParam = 0.1;
constexpr double kHelmholtzParam = 0.1;

struct SourceData {
  dashmm::Point position;
  double charge;
  int index;
};

struct TargetData {
  dashmm::Point position;
  std::complex<double> phi;
  int index;
};

// TODO eventually this should be exported to another header file so that
// anyone needing to read the file can get at this information.
struct FileHeader {
  int n_sources;
  int n_targets;
  int refinement_limit;
  double yukawa_param;
  double helmholtz_param;
  bool has_laplace;
  bool has_yukawa;
  bool has_helmholtz;
};

struct FileTargetData {
  dashmm::Point position;
  std::complex<double> phi;
  std::complex<double> phi_laplace;
  std::complex<double> phi_yukawa;
  std::complex<double> phi_helmholtz;
  int index;
};

// Here we create the evaluator objects that we shall need in this demo.
// These must be instantiated before the call to dashmm::init so that they
// might register the relevant actions with the runtime system.
//
// It is more typical to have only one or possibly two evaluators. The
// proliferation here is because this demo program allows for many use-cases.
dashmm::Evaluator<SourceData, TargetData,
                  dashmm::Laplace, dashmm::Direct> laplace_direct{};
dashmm::Evaluator<SourceData, TargetData,
                  dashmm::Yukawa, dashmm::Direct> yukawa_direct{};
dashmm::Evaluator<SourceData, TargetData,
                  dashmm::Helmholtz, dashmm::Direct> helmholtz_direct{};

// This type collects the input arguments to the program.
struct InputArguments {
  std::string type;
  std::string kernel;
};

// Print usage information.
void print_usage(char *progname) {
  fprintf(stdout, "Usage: %s [OPTIONS]\n\n"
          "Options available: [possible/values] (default value)\n"
          "--kernel=[laplace/yukawa/helmholtz]  method to use (fmm)\n"
          "--data=[cube/sphere]                 distribution type (cube)\n"
          , progname);
}

// Parse the command line arguments, overiding any defaults at the request of
// the user.
int read_arguments(int argc, char **argv, InputArguments &retval) {
  //Set defaults
  retval.type = std::string{"cube"};
  retval.kernel = std::string{"laplace"};

  int opt = 0;
  static struct option long_options[] = {
    {"data", required_argument, 0, 'w'},
    {"kernel", required_argument, 0, 'k'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
  };

  int long_index = 0;
  while ((opt = getopt_long(argc, argv, "w:k:h",
                            long_options, &long_index)) != -1) {
    std::string verifyarg{};
    switch (opt) {
    case 'w':
      retval.type = optarg;
      break;
    case 'k':
      retval.kernel = optarg;
      break;
    case 'h':
      print_usage(argv[0]);
      return -1;
    case '?':
      return -1;
    }
  }

  //test the inputs
  if (retval.type != "cube" && retval.type != "sphere"
        && retval.type != "plummer") {
    fprintf(stderr, "Usage ERROR: unknown source type '%s'\n",
            retval.type.c_str());
    return -1;
  }

  //print out summary
  fprintf(stdout, "Creating DASHMM test points:\n");
  fprintf(stdout, "%d sources and %d targets in a %s distribution\n",
          kNSources, kNTargets, retval.type.c_str());
  fprintf(stdout, "kernel: %s\n\n", retval.kernel.c_str());

  return 0;
}

// Pick the positions in a cube with a uniform distribution
dashmm::Point pick_cube_position() {
  double pos[3];
  pos[0] = (double)rand() / RAND_MAX - 0.5;
  pos[1] = (double)rand() / RAND_MAX - 0.5;
  pos[2] = (double)rand() / RAND_MAX - 0.5;
  return dashmm::Point{pos[0], pos[1], pos[2]};
}

// Pick a point from the surface of the sphere, with a uniform distribution
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

// Pick a random charge vale
double pick_charge() {
  double retval = (double)rand() / RAND_MAX + 1.0;
  if (rand() % 2) {
    retval *= -1.0;
  }
  return retval;
}

// Set the source data
void set_sources(SourceData *sources, int source_count,
                 std::string type) {
  if (type == std::string{"cube"}) {
    for (int i = 0; i < source_count; ++i) {
      sources[i].position = pick_cube_position();
      sources[i].charge = pick_charge();
      sources[i].index = i;
    }
  } else if (type == std::string{"sphere"}) {
    //Sphere
    for (int i = 0; i < source_count; ++i) {
      sources[i].position = pick_sphere_position();
      sources[i].charge = pick_charge();
      sources[i].index = i;
    }
  }
}

// Prepare the source data. This means both allocate the dashmm::Array,
// as well as creating the source data, and putting the source data into the
// global address space.
dashmm::Array<SourceData> prepare_sources(InputArguments &args) {
  SourceData *sources{nullptr};
  sources = new SourceData[kNSources];
  set_sources(sources, kNSources, args.type);

  dashmm::Array<SourceData> retval{};
  int err = retval.allocate(kNSources);
  assert(err == dashmm::kSuccess);
  err = retval.put(0, kNSources, sources);
  assert(err == dashmm::kSuccess);

  delete [] sources;

  return retval;
}

// Set the target data
void set_targets(TargetData *targets, int target_count,
                 std::string type) {
  if (type == std::string{"cube"}) {
    //Cube
    for (int i = 0; i < target_count; ++i) {
      targets[i].position = pick_cube_position();
      targets[i].phi = 0.0;
      targets[i].index = i;
    }
  } else if (type == std::string{"sphere"}) {
    //Sphere
    for (int i = 0; i < target_count; ++i) {
      targets[i].position = pick_sphere_position();
      targets[i].phi = 0.0;
      targets[i].index = i;
    }
  }
}

// Prepare the target data. This means both allocate the dashmm::Array,
// as well as creating the target data, and putting the target data into the
// global address space.
dashmm::Array<TargetData> prepare_targets(InputArguments &args) {
  TargetData *targets{nullptr};
  targets = new TargetData[kNTargets];
  set_targets(targets, kNTargets, args.type);

  dashmm::Array<TargetData> retval{};
  int err = retval.allocate(kNTargets);
  assert(err == dashmm::kSuccess);
  err = retval.put(0, kNTargets, targets);
  assert(err == dashmm::kSuccess);

  if (kNTargets) {
    delete [] targets;
  }

  return retval;
}


// Save the results to file
void save_to_file(dashmm::Array<SourceData> source_handle,
                  dashmm::Array<TargetData> target_handle,
                  const std::string &kernel,
                  const std::string &shape) {
  // Collect the data
  auto sources = source_handle.collect();
  auto targets = target_handle.collect();

  // Sort by the index
  SourceData *slocal = sources.get();
  std::sort(slocal, &slocal[kNSources],
            [](const SourceData &a, const SourceData &b) -> bool {
              return a.index < b.index;
            });
  TargetData *tlocal = targets.get();
  std::sort(tlocal, &tlocal[kNTargets],
            [](const TargetData &a, const TargetData &b) -> bool {
              return a.index < b.index;
            });

  // Open the file
  char fname[300];
  sprintf(fname, "prepared.%s.%s.dat", kernel.c_str(), shape.c_str());
  FILE *ofd = fopen(fname, "wb");
  assert(ofd != nullptr);

  // Prepare the header, and write it
  FileHeader head{
    kNSources,
    kNTargets,
    kRefinementLimit,
    kYukawaParam,
    kHelmholtzParam,
    false,
    false,
    false
  };
  std::function<void(const TargetData *, FileTargetData *)> saver{};
  if (kernel == std::string{"laplace"}) {
    head.has_laplace = true;
    saver = [](const TargetData *a, FileTargetData *b) {
      b->position = a->position;
      b->index = a->index;
      b->phi_laplace = a->phi;
    };
  } else if (kernel == std::string{"yukawa"}) {
    head.has_yukawa = true;
    saver = [](const TargetData *a, FileTargetData *b) {
      b->position = a->position;
      b->index = a->index;
      b->phi_yukawa = a->phi;
    };
  } else if (kernel == std::string{"helmholtz"}) {
    head.has_helmholtz = true;
    saver = [](const TargetData *a, FileTargetData *b) {
      b->position = a->position;
      b->index = a->index;
      b->phi_helmholtz = a->phi;
    };
  }
  assert(1 == fwrite(&head, sizeof(head), 1, ofd));

  // write the sources
  assert(kNSources == fwrite(slocal, sizeof(SourceData), kNSources, ofd));

  // write the targets
  for (int i = 0; i < kNSources; ++i) {
    FileTargetData loop{};
    saver(&targets[i], &loop);
    assert(1 == fwrite(&loop, sizeof(loop), 1, ofd));
  }

  // close the file
  fclose(ofd);
}


// The main driver routine that performes the test of evaluate()
void perform_evaluation_test(InputArguments args) {
  srand(123456);

  dashmm::Array<SourceData> source_handle = prepare_sources(args);
  dashmm::Array<TargetData> target_handle = prepare_targets(args);

  //Perform the evaluation
  int err{0};
  if (args.kernel == std::string{"laplace"}) {
    dashmm::Direct<SourceData, TargetData, dashmm::Laplace> direct{};
    err = laplace_direct.evaluate(source_handle, target_handle,
                                  kRefinementLimit, direct, 3,
                                  std::vector<double>{});
    assert(err == dashmm::kSuccess);
  } else if (args.kernel == std::string{"yukawa"}) {
    dashmm::Direct<SourceData, TargetData, dashmm::Yukawa> direct{};
    std::vector<double> kernelparms(1, kYukawaParam);
    err = yukawa_direct.evaluate(source_handle, target_handle,
                                 kRefinementLimit, direct, 3, kernelparms);
    assert(err == dashmm::kSuccess);
  } else if (args.kernel == std::string{"helmholtz"}) {
    dashmm::Direct<SourceData, TargetData, dashmm::Helmholtz> direct{};
    std::vector<double> kernelparms(1, kHelmholtzParam);
    err = helmholtz_direct.evaluate(source_handle, target_handle,
                                    kRefinementLimit, direct, 3, kernelparms);
    assert(err == dashmm::kSuccess);
  }

  // Now save the data to file
  save_to_file(source_handle, target_handle, args.kernel, args.type);

  //free up resources
  err = source_handle.destroy();
  assert(err == dashmm::kSuccess);
  err = target_handle.destroy();
  assert(err == dashmm::kSuccess);
}

// Program entrypoint
int main(int argc, char **argv) {
  auto err = dashmm::init(&argc, &argv);
  assert(err == dashmm::kSuccess);

  if (dashmm::get_num_ranks() > 1) {
    fprintf(stderr, "ERROR: Only designed for one rank.");
    return -1;
  }

  InputArguments inputargs;
  int usage_error = read_arguments(argc, argv, inputargs);

  if (!usage_error) {
    perform_evaluation_test(inputargs);
  }

  err = dashmm::finalize();
  assert(err == dashmm::kSuccess);

  return 0;
}
