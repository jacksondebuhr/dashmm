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


#ifndef __DASHMM_SPECIAL_FUNCTION_H__
#define __DASHMM_SPECIAL_FUNCTION_H__


/// \file
/// \brief Declaration of special functions


#include <cmath>
#include <cassert>
#include "dashmm/types.h"


namespace dashmm {

/// Offset of coefficient (*)_n^m, where 0 <= m <= n
inline int midx(const int n, const int m) {
  return n * (n + 1) / 2 + m;
}

/// Offset of coefficient (*)_n^m, where -n <= m <= n
inline int lidx(const int n, const int m) {
  return n * (n + 1) + m;
}

inline int didx(const int n, const int mp, const int m) {
  return n * (n + 1) * (4 * n - 1) / 6 + mp * (2 * n + 1) + n + m;
}

inline int hidx(const int n, const int mp, const int m) {
  return n * (2 * n - 1) * (2 * n + 1) / 3 + (mp + n) * (2 * n + 1) + n + m;
}

inline int sidx(int n, int m, int np, int p) {
  return midx(n, m) * (p + 1) + np;
}

inline double pow_m1(const int m) {
  return (m % 2 ? -1.0 : 1.0);
}

/// Compute Legendre polynomial P_n^m(x), where |x| <= 1, 0 <= m <= n
void legendre_Plm(int n, double x, double *P);

/// Compute scaled Legendre polynomial scale^n P_n^m(x) where |x| > 1
void legendre_Plm_gt1_scaled(int nb, double x, double scale, double *P);

/// Compute scaled Legendre polynomial scale^n R_n^{|m|}(x), where
/// R_n^{|m|}(x) i^{n+|m|} = P_n^{|m|}(ix)
void legendre_Plm_evan_scaled(int nb, double x, double scale, double *P);

/// Compute scaled Legendre polynomial (-i * scale)^n P_n^m(x)
void legendre_Plm_prop_scaled(int nb, double x, double scale, dcomplex_t *P);

/// Gamma function
double Gamma(double x);

/// Cylindrical Bessel function of the first kind
int bessel_Jn(int nb, double alpha, double x, double *B);

/// Cylindrical Bessel function of the second kind
int bessel_Yn(int nb, double alpha, double x, double *B);

/// Modified cylindrical Bessel function of the first kind
int bessel_In(int nb, double alpha, double x, int ize, double *B);

/// Spherical Bessel function of the first kind
void bessel_jn_scaled(int nb, double x, double scale, double *B);

/// Modified spherical Bessel function of the first kind
void bessel_in_scaled(int nb, double x, double scale, double *B);

/// Modified spherical Bessel function of the second kind
void bessel_kn_scaled(int nb, double x, double scale, double *B);

/// Spherical Hankel function of the first kind
void bessel_hn_scaled(int nb, double x, double scale, dcomplex_t *B);


} // namespace dashmm


#endif // __DASHMM_SPECIAL_FUNCTION_H__
