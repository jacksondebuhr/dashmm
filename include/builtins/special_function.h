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

#ifndef __DASHMM_SPECIAL_FUNCTION_H__
#define __DASHMM_SPECIAL_FUNCTION_H__

/// \file include/builtins/special_function.h
/// \brief Declaration of special functions 

#include <cmath>
#include <cassert>

namespace dashmm {

inline int midx(const int n, const int m) {
  return n * (n + 1) / 2 + m;
}

inline int didx(const int n, const int mp, const int m) {
  return n * (n + 1) * (4 * n - 1) / 6 + mp * (2 * n + 1) + n + m;
}

inline int sidx(int n, int m, int np, int p) {
  return midx(n, m) * (p + 1) + np;
}

inline double pow_m1(const int m) {
  return (m % 2 ? -1.0 : 1.0);
}

void legendre_Plm(int n, double x, double *P);

void legendre_Plm_gt1_scaled(int nb, double x, double scale, double *P); 

double Gamma(double x); 

int bessel_In(int nb, double alpha, double x, int ize, double *B); 

void bessel_in_scaled(int nb, double x, double scale, double *B); 

void bessel_kn_scaled(int nb, double x, double scale, double *B); 

} // namespace dashmm 

#endif // __DASHMM_SPECIAL_FUNCTION_H__
