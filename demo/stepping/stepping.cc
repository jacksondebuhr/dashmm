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

#include "dashmm/dashmm.h"


// NOTE: We are working in units where G = 1.
struct Particle {
  dashmm::Point position;
  double charge;
  double velocity[3];
  double acceleration[3];
  int index;
};


struct InputArguments {
  int count;
  int refinement_limit;
  int steps;
  std::string output;
};


void update_particles(Particle *P, const size_t count, const size_t offset,
                      const double *dt);


// Here we create the evaluator objects that we shall need in this demo.
// This must be instantiated before the call to dashmm::init so that they
// might register the relevant actions with the runtime system.
dashmm::Evaluator<Particle, Particle,
                  dashmm::LaplaceCOMAcc, dashmm::BH> bheval{ };

// We also create the ArrayMapAction object before dashmm::init so it too
// can register the relevant actions with the runtime system.
dashmm::ArrayMapAction<Particle, double> update_action{update_particles};


void print_usage(char *progname) {
  fprintf(stdout, "Usage: %s [OPTIONS]\n\n"
"Options available: [possible/values] (default value)\n"
"  --nsources=num               number of source points to generate (10000)\n"
"  --threshold=num              source and target tree partition refinement\n"
"                                 limit (40)\n"
"  --nsteps=num                 number of steps to take (100)\n"
"  --output=file                specify file for output (disabled)\n",
          progname);
}


int read_arguments(int argc, char **argv, InputArguments &retval) {
  // Set defaults
  retval.count = 10000;
  retval.refinement_limit = 40;
  retval.steps = 100;
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
  if (hpx_get_my_rank() == 0) {
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


inline double getticks(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double) (tv.tv_sec * 1e6 + tv.tv_usec);
}

inline double elapsed(double t1, double t0) {
  return (double) (t1 - t0);
}


dashmm::Point pick_sphere_position() {
  double r = (double)rand() / RAND_MAX;
  double ctheta = 2.0 * (double)rand() / RAND_MAX - 1.0;
  double stheta = sqrt(1.0 - ctheta * ctheta);
  double phi = 2.0 * 3.1415926535 * (double)rand() / RAND_MAX;
  double pos[3];
  pos[0] = r * stheta * cos(phi);
  pos[1] = r * stheta * sin(phi);
  pos[2] = r * ctheta;
  return dashmm::Point{pos[0], pos[1], pos[2]};
}


double pick_mass() {
  double retval = (double)rand() / RAND_MAX + 1.0;
  return retval;
}


double set_sources(Particle *sources, int count) {
  double m_tot{0.0};
  for (int i = 0; i < count; ++i) {
    sources[i].position = pick_sphere_position();
    sources[i].charge = pick_mass();
    m_tot += sources[i].charge;
    sources[i].velocity[0] = 0.0;
    sources[i].velocity[1] = 0.0;
    sources[i].velocity[2] = 0.0;
    sources[i].acceleration[0] = 0.0;
    sources[i].acceleration[1] = 0.0;
    sources[i].acceleration[2] = 0.0;
    sources[i].index = i;
  }
  return m_tot;
}


void output_results(const std::string &fname, const Particle *sources,
                    int count) {
  if (hpx_get_my_rank()) return;

  FILE *ofd = fopen(fname.c_str(), "w");
  assert(ofd != nullptr);

  for (int i = 0; i < count; ++i) {
    fprintf(ofd, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %d\n",
            sources[i].position[0], sources[i].position[1],
            sources[i].position[2], sources[i].velocity[0],
            sources[i].velocity[1], sources[i].velocity[2],
            sources[i].acceleration[0], sources[i].acceleration[1],
            sources[i].acceleration[2], sources[i].index);
  }

  fclose(ofd);
}


void update_particles(Particle *P, const size_t count, const size_t offset,
                      const double *dt) {
  for (size_t i = 0; i < count; ++i) {
    double x[3] = {P[i].position[0], P[i].position[1], P[i].position[2]};
    for (int j = 0; j < 3; ++j) {
      // NOTE: the minus sign is because the output of LaplaceCOMAcc is for
      // a positive potential. So we add the minus sign here.
      // NOTE: we work in units where G = 1.
      // NOTE: the point of this demo is not the time integrator, hence the
      // simplistic update.
      x[j] += P[i].velocity[j] * (*dt)
              - 0.5 * P[i].acceleration[j] * (*dt) * (*dt);
      P[i].velocity[j] -= P[i].acceleration[j] * (*dt);
    }
    P[i].position = dashmm::Point(x);
  }
}


void perform_time_stepping(InputArguments args) {
  srand(123456);

  // create some arrays
  Particle *sources{nullptr};
  double m_tot{1.0}; // This is given a value to avoid divbyzero below
  if (args.count) {
    sources = reinterpret_cast<Particle *>(
      new char [sizeof(Particle) * args.count]);
    m_tot = set_sources(sources, args.count);
  }

  // prep source Array
  dashmm::Array<Particle> source_handle{ };
  int err = source_handle.allocate(args.count);
  assert(err == dashmm::kSuccess);
  err = source_handle.put(0, args.count, sources);
  assert(err == dashmm::kSuccess);

  // Compute a reasonable dt:
  // This is chosen so that one dynamical time is roughly 200 steps.
  // NOTE: this will only be correct on rank 0, but it will also only be used
  // by rank zero.
  double dt{sqrt(1.0 / m_tot) / 200.0};

  // Clear timing information
  double t_eval{0.0};
  double t_update{0.0};

  // Prototypes for the expansion and method
  dashmm::LaplaceCOMAcc<Particle, Particle> expansion{
    dashmm::Point{0.0, 0.0, 0.0}, 0, dashmm::kNoRoleNeeded
  };
  dashmm::BH<Particle, Particle, dashmm::LaplaceCOMAcc> method{0.6};

  // Time-stepping
  for (int step = 0; step < args.steps; ++step) {
    double t0 = getticks();
    err = bheval.evaluate(source_handle, source_handle, args.refinement_limit,
                          method, expansion);
    assert(err == dashmm::kSuccess);
    double t1 = getticks();

    // Now update the positions based on the velocity
    source_handle.map(update_action, &dt);
    double t2 = getticks();

    // Collect timing
    t_eval += elapsed(t1, t0);
    t_update += elapsed(t2, t1);
  }

  // Report on loop
  fprintf(stdout, "Evaluation took %lg [us]\n", t_eval);
  fprintf(stdout, "Update took %lg [us]\n", t_update);

  // Output if the user has selected this option
  if (!args.output.empty()) {
    err = source_handle.get(0, args.count, sources);
    assert(err == dashmm::kSuccess);
    output_results(args.output, sources, args.count);
  }

  // free up resources
  err = source_handle.destroy();
  assert(err == dashmm::kSuccess);

  delete [] sources;
}


int main(int argc, char **argv) {
  auto err = dashmm::init(&argc, &argv);
  assert(err == dashmm::kSuccess);


  InputArguments inputargs;
  int usage_error = read_arguments(argc, argv, inputargs);


  if (!usage_error) {
    perform_time_stepping(inputargs);
  }


  err = dashmm::finalize();
  assert(err == dashmm::kSuccess);

  return 0;
}
