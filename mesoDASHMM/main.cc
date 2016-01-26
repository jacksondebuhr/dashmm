#include <algorithm>
#include <chrono>
#include <cstdio>
#include <functional>
#include <random>
#include <string>
#include <vector>
#include <getopt.h>

#include "dashmm.h"

using namespace std::chrono;


/// randomly generates particles using number given by the user
void data_generator(std::vector<dashmm::Source> &sources,
                    std::vector<dashmm::Target> &targets,
                    std::string datatype, 
                    std::string test_case) {
#if 0
  std::mt19937 engine(1);
  std::uniform_real_distribution<double> m_distrib(-1.0, 1.0);
  auto m_roll = std::bind(m_distrib, engine);

  if (datatype.compare("cube") == 0) {
    std::uniform_real_distribution<double> p_distrib(0.0, 1.0);
    auto p_roll = std::bind(p_distrib, engine);

    for (auto i = sources.begin(); i != sources.end(); ++i) {
      (*i) = dashmm::Source{p_roll(), p_roll(), p_roll(), m_roll()};
    }

    for (auto i = targets.begin(); i != targets.end(); ++i) {
      (*i) = dashmm::Target{p_roll(), p_roll(), p_roll()};
    }
  } else if (datatype.compare("sphere") == 0) {
    //NOTE: this is not uniform over the sphere. There will be more points
    // near the poles. If uniform was the intent, this should be drawn from
    // cos(theta) between -1 and 1, and not theta itself.
    std::uniform_real_distribution<double> theta_distrib(0.0, M_PI);
    std::uniform_real_distribution<double> phi_distrib(0.0, 2 * M_PI);
    auto theta_roll = std::bind(theta_distrib, engine);
    auto phi_roll = std::bind(phi_distrib, engine);

    for (auto i = sources.begin(); i != sources.end(); ++i) {
      double theta = theta_roll();
      double phi = phi_roll();
      (*i) = dashmm::Source{sin(theta) * cos(phi),
                            sin(theta) * sin(phi),
                            cos(theta), m_roll()};
    }

    for (auto i = targets.begin(); i != targets.end(); ++i) {
      double theta = theta_roll();
      double phi = phi_roll();
      (*i) = dashmm::Target{sin(theta) * cos(phi),
                            sin(theta) * sin(phi),
                            cos(theta)};
    }
  }
#endif 
  srand(123456); 
  bool use_negative = test_case == std::string{"fmm"}; 
  if (datatype.compare("cube") == 0) {
    for (auto i = sources.begin(); i != sources.end(); ++i) {
      double x = (double) rand() / RAND_MAX; 
      double y = (double) rand() / RAND_MAX; 
      double z = (double) rand() / RAND_MAX; 
      double m = (double) rand() / RAND_MAX + 1.0; 
      if (use_negative && (rand() % 2)) {
        m *= -1;
      }
      (*i) = dashmm::Source{x, y, z, m}; 
    }

    for (auto i = targets.begin(); i != targets.end(); ++i) {
      double x = (double) rand() / RAND_MAX; 
      double y = (double) rand() / RAND_MAX; 
      double z = (double) rand() / RAND_MAX; 
      (*i) = dashmm::Target{x, y, z}; 
    }    
  } else if (datatype.compare("sphere") == 0) {
    for (auto i = sources.begin(); i != sources.end(); ++i) {
      double r = 1.0; 
      double ctheta = 2.0 * (double)rand() / RAND_MAX - 1.0; 
      double stheta = sqrt(1.0 - ctheta * ctheta); 
      double phi = 2.0 * 3.1415926535 * (double) rand() / RAND_MAX; 
      double x = r * stheta * cos(phi); 
      double y = r * stheta * sin(phi); 
      double z = r * ctheta; 
      double m = (double) rand() / RAND_MAX + 1.0; 
      if (use_negative && (rand() % 2)) {
        m *= -1;
      }
      (*i) = dashmm::Source{x, y, z, m}; 
    }

    for (auto i = targets.begin(); i != targets.end(); ++i) {
      double r = 1.0; 
      double ctheta = 2.0 * (double)rand() / RAND_MAX - 1.0; 
      double stheta = sqrt(1.0 - ctheta * ctheta); 
      double phi = 2.0 * 3.1415926535 * (double) rand() / RAND_MAX; 
      double x = r * stheta * cos(phi); 
      double y = r * stheta * sin(phi); 
      double z = r * ctheta; 
      (*i) = dashmm::Target{x, y, z}; 
    }    
  } else {
    // Plummer 
    int count = sources.size(); 
    for (auto i = sources.begin(); i != sources.end(); ++i) {
      double unif = (double)rand() / RAND_MAX; 
      double r = 1.0 / sqrt(pow(unif, -2.0 / 3.0) - 1); 
      double ctheta = 2.0 * (double) rand() / RAND_MAX - 1.0; 
      double stheta = sqrt(1.0 - ctheta * ctheta); 
      double phi = 2.0 * 3.1415926535 * (double) rand() / RAND_MAX; 
      double x = r * stheta * cos(phi); 
      double y = r * stheta * sin(phi); 
      double z = r * ctheta; 
      double m = 100.0 / count; 
      (*i) = dashmm::Source{x, y, z, m}; 
    }

    for (auto i = targets.begin(); i != targets.end(); ++i) {
      double unif = (double)rand() / RAND_MAX; 
      double r = 1.0 / sqrt(pow(unif, -2.0 / 3.0) - 1); 
      double ctheta = 2.0 * (double) rand() / RAND_MAX - 1.0; 
      double stheta = sqrt(1.0 - ctheta * ctheta); 
      double phi = 2.0 * 3.1415926535 * (double) rand() / RAND_MAX; 
      double x = r * stheta * cos(phi); 
      double y = r * stheta * sin(phi); 
      double z = r * ctheta; 
      (*i) = dashmm::Target{x, y, z}; 
    }    
  }
}

/// calculates the potential using the direct sum method
void Particle_direct(std::vector<dashmm::Source> &sources,
                     std::vector<dashmm::Target> &targets,
                     std::vector<double> &relerr,
                     std::vector<double> &exacts) {
  std::mt19937 engine(
      std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> uniformdistrib(0.0, 1.0);
  auto useroll = std::bind(uniformdistrib, engine);

  double odds = 400.0 / targets.size();

  relerr.clear();
  exacts.clear();

  for (size_t i = 0; i < targets.size(); ++i) {
    double testit = useroll();
    if (testit > odds) continue;

    double phi = 0;
    dashmm::Point tloc = targets[i].position();
    for (size_t j = 0; j < sources.size(); ++j) {
      dashmm::Point sloc = sources[j].position();
      double distance = (sloc - tloc).norm();
      if (distance > 0)
        phi += sources[j].charge() / (distance);
    }
    targets[i].set_direct(phi);

    relerr.push_back(targets[i].direct().real() - targets[i].phi().real());
    exacts.push_back(targets[i].direct().real());
  }
}


/// gives the names of the input parameters for main
void print_usage(std::string progname) {
  fprintf(stdout, "Usage: %s --method=[fmm/bh] "
          "--nsources=num --ntargets=num "
          "--threshold=num --data=[cube/sphere] --verify=[yes/no]\n",
          progname.c_str());
}

int main(int argc, char **argv){
  // default values
  std::string method{"fmm"};
  int nsources = 1000;
  int ntargets = 1000;
  int threshold = 40;
  std::string datatype{"cube"};
  std::string verify{"yes"};

  int opt = 0;
  static struct option long_options[] = {
    {"method", required_argument, 0, 'm'},
    {"nsources", required_argument, 0, 's'},
    {"ntargets", required_argument, 0, 't'},
    {"threshold", required_argument, 0, 'l'},
    {"data", required_argument, 0, 'd'},
    {"verify", required_argument, 0, 'v'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
  };

  int long_index = 0;
  while ((opt = getopt_long(argc, argv, "m:s:t:l:d:v:h",
                            long_options, &long_index)) != -1) {
    switch (opt) {
    case 'm':
      method = optarg;
      break;
    case 's':
      nsources = atoi(optarg);
      break;
    case 't':
      ntargets = atoi(optarg);
      break;
    case 'l':
      threshold = atoi(optarg);
      break;
    case 'd':
      datatype = optarg;
      break;
    case 'v':
      verify = optarg;
      break;
    case 'h':
      print_usage(std::string{argv[0]});
      return -1;
    }
  }

  dashmm::init();
  std::vector<dashmm::Source> sources(nsources);
  std::vector<dashmm::Target> targets(ntargets);
  data_generator(sources, targets, datatype, method);

  dashmm::Expansion *exp_select{nullptr};
  dashmm::Method *method_select{nullptr};
  if (method.compare("fmm") == 0) {
    exp_select = dashmm::fmm_expansion(9);
    method_select = dashmm::fmm_method();
  } else {
    exp_select = dashmm::bh_expansion();
    method_select = dashmm::bh_method(0.65);
  }

  auto t1 = high_resolution_clock::now();
  dashmm::evaluate(sources, targets, threshold, method_select, exp_select);
  auto t2 = high_resolution_clock::now();


  duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
  fprintf(stdout, "%d %d %d %lg\n",
          nsources, ntargets, threshold, time_span.count());

  if (verify.compare("yes") == 0) {
    std::vector<double> relerr{};
    std::vector<double> exacts{};
    Particle_direct(sources, targets, relerr, exacts);

    //compute the standard error expression
    double numerator{0.0};
    double denominator{0.0};
    for (size_t i = 0; i < relerr.size(); ++i) {
      numerator += relerr[i] * relerr[i];
      denominator += exacts[i] * exacts[i];
    }
    double errorestimate = sqrt(numerator / denominator);

    //do this last, as we are scrambling the input
    double maxabserr{0.0};
    size_t index_of_max{0};
    for (size_t i = 0; i < relerr.size(); ++i) {
      if (fabs(relerr[i]) > maxabserr) {
        maxabserr = fabs(relerr[i]);
        index_of_max = i;
      }
      //fprintf(stdout, "%lg\n", fabs(relerr[i]) / fabs(exacts[i]));
    }
    double maxrelerr = fabs(relerr[index_of_max]) / fabs(exacts[index_of_max]);

    fprintf(stdout, "Number of points tested: %lu\n", relerr.size());
    fprintf(stdout, "Error characterization: %lg (max relative error %lg)\n",
                    errorestimate, maxrelerr);
  }


  dashmm::finalize();
  return 0;
}
