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
#include <map>
#include <memory>
#include <string>
#include "dashmm/dashmm.h"

#include "../common/common.h"

// Alias to shorten name somewhat
using Source = FileSourceData;
using Target = FileTargetData;

constexpr double kBHTheta = 0.6;

// Here we create the evaluator objects that we shall need in this demo.
// These must be instantiated before the call to dashmm::init so that they
// might register the relevant actions with the runtime system.
//
// It is more typical to have only one or possibly two evaluators. The
// proliferation here is because this demo program allows for many use-cases.
dashmm::Evaluator<Source, Target,
                  dashmm::LaplaceCOM, dashmm::BH> laplace_bh{};
dashmm::Evaluator<Source, Target,
                  dashmm::Laplace, dashmm::FMM> laplace_fmm{};
dashmm::Evaluator<Source, Target,
                  dashmm::Laplace, dashmm::FMM97> laplace_fmm97{};
dashmm::Evaluator<Source, Target,
                  dashmm::Yukawa, dashmm::FMM97> yukawa_fmm97{};
dashmm::Evaluator<Source, Target,
                  dashmm::Helmholtz, dashmm::FMM97> helmholtz_fmm97{};

// This type collects the input arguments to the program.
struct InputArguments {
  std::string datafile;
  int refinement_limit;
  std::string method;
  std::string kernel;
  int accuracy;
};

// Print usage information.
void print_usage(char *progname) {
  fprintf(stdout, "Usage: %s [OPTIONS]\n\n"
          "Options available: [possible/values] (default value)\n"
          "--method=[fmm/fmm97/bh]     method to use (fmm)\n"
          "--data=[filename]\n         source/target datafile\n"
          "--threshold=num             "
          "source and target tree partition refinement limit (40)\n"
          "--accuracy=num              "
          "number of digits of accuracy for fmm (3)\n"
          "--kernel=[laplace/yukawa/helmholtz]   "
          "particle interaction type (laplace)\n"
          , progname);
}

// Parse the command line arguments, overiding any defaults at the request of
// the user.
int read_arguments(int argc, char **argv, InputArguments &retval) {
  //Set defaults
  retval.datafile = std::string{""};
  retval.refinement_limit = 40;
  retval.method = std::string{"fmm97"};
  retval.kernel = std::string{"laplace"};
  retval.accuracy = 3;

  int opt = 0;
  static struct option long_options[] = {
    {"method", required_argument, 0, 'm'},
    {"data", required_argument, 0, 'w'},
    {"threshold", required_argument, 0, 'l'},
    {"accuracy", required_argument, 0, 'a'},
    {"kernel", required_argument, 0, 'k'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
  };

  int long_index = 0;
  while ((opt = getopt_long(argc, argv, "m:s:w:t:g:l:v:a:k:h",
                            long_options, &long_index)) != -1) {
    std::string verifyarg{};
    switch (opt) {
    case 'm':
      retval.method = optarg;
      break;
    case 'w':
      retval.datafile = optarg;
      break;
    case 'l':
      retval.refinement_limit = atoi(optarg);
      break;
    case 'a':
      retval.accuracy = atoi(optarg);
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
  if (retval.datafile == "") {
    fprintf(stderr, "Usage ERROR: must specify input data file\n");
    return -1;
  }


  if (retval.method != "bh" && retval.method != "fmm"
      && retval.method != "fmm97") {
    fprintf(stderr, "Usage ERROR: unknown method '%s'\n",
            retval.method.c_str());
    return -1;
  }

  if (retval.kernel == "laplace" && retval.method == "fmm97") {
    if (retval.accuracy != 3 && retval.accuracy != 6) {
      fprintf(stderr, "Usage ERROR: only 3-/6-digit accuracy supported"
              " for laplace kernel using fmm97\n");
      return -1;
    }
  }

  if (retval.kernel == "yukawa") {
    if (retval.method != "fmm97") {
      fprintf(stderr, "Usage ERROR: yukawa kernel must use fmm97\n");
      return -1;
    } else if (retval.accuracy != 3 && retval.accuracy != 6) {
      fprintf(stderr, "Usage ERROR: only 3-/6-digit accuracy supported"
              " for yukawa kernel using fmm97\n");
      return -1;
    }
  }

  if (retval.kernel == "helmholtz") {
    if (retval.method != "fmm97") {
      fprintf(stderr, "Usage ERROR: helmholtz kernel must use fmm97\n");
      return -1;
    } else if (retval.accuracy != 3 && retval.accuracy != 6) {
      fprintf(stderr, "Usage ERROR: only 3-/6-digit accuracy supported"
              " for helmholtz kernel using fmm97\n");
      return -1;
    }
  }

  //print out summary
  if (dashmm::get_my_rank() == 0) {
    fprintf(stdout, "Testing DASHMM:\n");
    fprintf(stdout, "using data from: '%s'\n", retval.datafile.c_str());
    fprintf(stdout, "method: %s \nthreshold: %d\nkernel: %s\n\n",
            retval.method.c_str(), retval.refinement_limit,
            retval.kernel.c_str());
  }

  return 0;
}

// Used to time the execution
inline double getticks(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double) (tv.tv_sec * 1e6 + tv.tv_usec);
}

inline double elapsed(double t1, double t0) {
  return (double) (t1 - t0);
}


FileHeader read_input_data(const std::string &fname,
                           Source **sources, Target **targets) {
  // read in the file
  if (hpx_get_my_rank() == 0) {
    FILE *ifd = fopen(fname.c_str(), "rb");
    assert(ifd != nullptr);

    FileHeader retval{};
    assert(1 == fread(&retval, sizeof(retval), 1, ifd));

    *sources = new Source[retval.n_sources];
    *targets = new Target[retval.n_targets];

    assert(retval.n_sources == (int)fread(*sources, sizeof(Source),
                                     retval.n_sources, ifd));
    assert(retval.n_targets == (int)fread(*targets, sizeof(Target),
                                     retval.n_targets, ifd));

    fclose(ifd);

    return retval;
  } else {
    return FileHeader{0, 0, 0, false, false, false, 1.0, 1.0};
  }
}


dashmm::Array<Source> prepare_sources(int n_sources, Source *sources) {
  dashmm::Array<Source> retval{};
  int err = retval.allocate(n_sources);
  assert(err == dashmm::kSuccess);
  err = retval.put(0, n_sources, sources);
  assert(err == dashmm::kSuccess);
  return retval;
}


dashmm::Array<Target> prepare_targets(int n_targets, Target *targets) {
  dashmm::Array<Target> retval{};
  int err = retval.allocate(n_targets);
  assert(err == dashmm::kSuccess);
  err = retval.put(0, n_targets, targets);
  assert(err == dashmm::kSuccess);
  return retval;
}


void compare_results(Target *targets, int target_count,
                     const std::string &kernel) {
  if (dashmm::get_my_rank()) return;

  std::function<double(const Target &)> normit{nullptr};
  std::function<double(const Target &)> denom{nullptr};
  std::function<double(const Target &)> maxrelerr{nullptr};
  if (kernel == "laplace") {
    normit = [](const Target &a) -> double {
      return std::norm(a.phi - a.phi_laplace);
    };
    denom = [](const Target &a) -> double {
      return std::norm(a.phi_laplace);
    };
    maxrelerr = [](const Target &a) -> double {
      return fabs((a.phi.real() - a.phi_laplace.real()) / a.phi_laplace.real());
    };
  } else if (kernel == "yukawa") {
    normit = [](const Target &a) -> double {
      return std::norm(a.phi - a.phi_yukawa);
    };
    denom = [](const Target &a) -> double {
      return std::norm(a.phi_yukawa);
    };
    maxrelerr = [](const Target &a) -> double {
      return fabs((a.phi.real() - a.phi_yukawa.real()) / a.phi_yukawa.real());
    };
  } else { //must be helmholtz
    normit = [](const Target &a) -> double {
      return std::norm(a.phi - a.phi_helmholtz);
    };
    denom = [](const Target &a) -> double {
      return std::norm(a.phi_helmholtz);
    };
    maxrelerr = [](const Target &a) -> double {
      return fabs((a.phi.real() - a.phi_helmholtz.real())
                        / a.phi_helmholtz.real());
    };
  }

  //now we loop over the exact results and compare
  double numerator = 0.0;
  double denominator = 0.0;
  double maxrel = 0.0;
  for (int i = 0; i < target_count; ++i) {
    numerator += normit(targets[i]);;
    denominator += denom(targets[i]);
    double localrel = maxrelerr(targets[i]);
    if (localrel > maxrel) {
      maxrel = localrel;
    }
  }
  fprintf(stdout, "Error for %d test points: %4.3e (max %4.3e)\n",
                  target_count, sqrt(numerator / denominator), maxrel);
}

// The main driver routine that performes the test of evaluate()
void perform_evaluation_test(InputArguments args) {
  srand(123456);

  // Read in the data
  Source *sources{nullptr};
  Target *targets{nullptr};
  FileHeader header = read_input_data(args.datafile, &sources, &targets);
  dashmm::Array<Source> source_handle = prepare_sources(header.n_sources,
                                                        sources);
  dashmm::Array<Target> target_handle = prepare_targets(header.n_targets,
                                                        targets);
  delete [] sources;
  delete [] targets;

  //Perform the evaluation
  double t0{};
  double tf{};
  int err{0};

  if (args.kernel == std::string{"laplace"}) {
    if (args.method == std::string{"bh"}) {
      dashmm::BH<Source, Target, dashmm::LaplaceCOM> method{kBHTheta};

      t0 = getticks();
      err = laplace_bh.evaluate(source_handle, target_handle,
                                args.refinement_limit, method,
                                args.accuracy, std::vector<double>{});
      assert(err == dashmm::kSuccess);
      tf = getticks();
    } else if (args.method == std::string{"fmm"}) {
      dashmm::FMM<Source, Target, dashmm::Laplace> method{};

      t0 = getticks();
      err = laplace_fmm.evaluate(source_handle, target_handle,
                                 args.refinement_limit, method,
                                 args.accuracy, std::vector<double>{});
      assert(err == dashmm::kSuccess);
      tf = getticks();
    } else if (args.method == std::string{"fmm97"}) {
      dashmm::FMM97<Source, Target, dashmm::Laplace> method{};

      t0 = getticks();
      err = laplace_fmm97.evaluate(source_handle, target_handle,
                                   args.refinement_limit, method,
                                   args.accuracy, std::vector<double>{});
      assert(err == dashmm::kSuccess);
      tf = getticks();
    }
  } else if (args.kernel == std::string{"yukawa"}) {
    if (args.method == std::string{"fmm97"}) {
      dashmm::FMM97<Source, Target, dashmm::Yukawa> method{};
      std::vector<double> kernelparms(1, header.yukawa_param);

      t0 = getticks();
      err = yukawa_fmm97.evaluate(source_handle, target_handle,
                                  args.refinement_limit, method,
                                  args.accuracy, kernelparms);
      assert(err == dashmm::kSuccess);
      tf = getticks();
    }
  } else if (args.kernel == std::string{"helmholtz"}) {
    if (args.method == std::string{"fmm97"}) {
      dashmm::FMM97<Source, Target, dashmm::Helmholtz> method{};
      std::vector<double> kernelparms(1, header.helmholtz_param);

      t0 = getticks();
      err = helmholtz_fmm97.evaluate(source_handle, target_handle,
                                     args.refinement_limit, method,
                                     args.accuracy, kernelparms);
      assert(err == dashmm::kSuccess);
      tf = getticks();
    }
  }

  fprintf(stdout, "Evaluation took %lg [us]\n", elapsed(tf, t0));

  // collect targets
  auto results = target_handle.collect();

  // run comparison
  compare_results(results.get(), header.n_targets, args.kernel);

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

  InputArguments inputargs;
  int usage_error = read_arguments(argc, argv, inputargs);

  if (!usage_error) {
    perform_evaluation_test(inputargs);
  }

  err = dashmm::finalize();
  assert(err == dashmm::kSuccess);

  return 0;
}
