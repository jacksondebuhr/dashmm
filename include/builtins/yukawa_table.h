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

#ifndef __DASHMM_YUKAWA_TABLE_H__
#define __DASHMM_YUKAWA_TABLE_H__


/// \file include/builtins/yukawa_table.h
/// \brief Declaration of precomputed tables for Yukawa


#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <vector>
#include <cassert>
#include "dashmm/types.h"

namespace dashmm {

struct yukawa_cmp {
  bool operator()(const double &a, const double &b) const {
    // The smallest gap between key values in the Yukawa rotation matrix map is
    // 0.01. The operator compares the first 6 digits and that should be
    // enough.
    double aa = floor(a * 1000000.0) / 100000.0;
    double bb = floor(b * 1000000.0) / 100000.0;
    if (aa == bb)
      return false;
    return aa < bb;
  }
};

using yukawa_map_t = std::map<double, double *, yukawa_cmp>;

class YukawaTable {
 public:
  YukawaTable(int n_digits, double size, double lambda); 
  ~YukawaTable(); 
  static const int maxlev; 
  int p() const {return p_;}
  int s() const {return s_;}
  double scale(int lev) const {return scale_ / pow(2, lev);}
  double lambda() const {return lambda_;} 
  const double *sqf() const {return sqf_;}
  const double *dmat_plus(double v) const {return dmat_plus_->at(v);}
  const double *dmat_minus(double v) const {return dmat_minus_->at(v);} 
  const double *m2m(double scale) const {
    return &m2m_[(p_ + 1) * (p_ + 1) * (p_ + 2) / 2 * level(scale)];
  }
  const double *l2l(double scale) const {
    return &l2l_[(p_ + 1) * (p_ + 1) * (p_ + 2) / 2 * level(scale)]; 
  }
  const double *x() const {return x_;} 
  const double *w() const {return w_;} 
  const int *m(double scale) const {return &m_[s_ * level(scale)];}
  const int *sm(double scale) const {return &sm_[(s_ + 1) * level(scale)];}
  int nexp(double scale) const {return nexp_[level(scale)];} 
  const int *f() const {return f_;}
  const int *smf(double scale) const {return &smf_[(s_ + 1) * level(scale)];}
  const dcomplex_t *ealphaj(double scale) const {return ealphaj_[level(scale)];}
  const dcomplex_t *xs(double scale) const {return xs_[level(scale)];}
  const dcomplex_t *ys(double scale) const {return ys_[level(scale)];}
  const double *zs(double scale) const {return zs_[level(scale)];}
  double size(double scale) const {return size_ * scale / scale_;}

private: 
  int p_; 
  int s_; 
  double lambda_; 
  double size_; 
  double scale_; // scaling factor of level 0 to avoid under-/over-flow 
  double *sqf_; 
  yukawa_map_t *dmat_plus_; 
  yukawa_map_t *dmat_minus_; 
  double *m2m_; 
  double *l2l_; 
  double *x_; 
  double *w_; 
  int *m_; 
  int *sm_; 
  int *nexp_; 
  int *f_; 
  int *smf_; 
  dcomplex_t **ealphaj_; 
  dcomplex_t **xs_; 
  dcomplex_t **ys_; 
  double **zs_; 

  int level(double scale) const {return log2(scale_ / scale);} 
  void generate_sqf(); 
  void generate_scaled_wigner_dmat(); 
  void generate_scaled_dmat_of_beta(double beta, double *dp, double *dm); 
  void generate_m2m(); 
  void generate_l2l(); 
}; 

using uYukawaTable = std::unique_ptr<YukawaTable>; 
extern uYukawaTable builtin_yukawa_table_; 

void legendre_Plm(int n, double x, double *P);
void legendre_Plm_gt1_scaled(int nb, double x, double scale, double *P); 
double Gamma(double x); 
int bessel_In(int nb, double alpha, double x, int ize, double *B); 
void bessel_in_scaled(int nb, double x, double scale, double *B); 
void bessel_kn_scaled(int nb, double x, double scale, double *B); 

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

void get_or_add_yukawa_table(int n_digits, double size, double lambda); 

} // namespace dashmm


#endif // __DASHMM_YUKAWA_TABLE_H__
