#ifndef __DASHMM_TEST_COMMON_COMMON_H__
#define __DASHMM_TEST_COMMON_COMMON_H__


#include "dashmm/dashmm.h"


struct FileHeader {
  int n_sources;
  int n_targets;
  int refinement_limit;
  bool has_laplace;
  bool has_yukawa;
  bool has_helmholtz;
  double yukawa_param;
  double helmholtz_param;
};

struct FileSourceData {
  dashmm::Point position;
  double charge;
  int index;
};

struct FileTargetData {
  dashmm::Point position;
  std::complex<double> phi;
  std::complex<double> phi_laplace;
  std::complex<double> phi_yukawa;
  std::complex<double> phi_helmholtz;
  int index;
};


#endif // __DASHMM_TEST_COMMON_COMMON_H__