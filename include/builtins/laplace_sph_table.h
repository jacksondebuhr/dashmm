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


#ifndef __DASHMM_LAPLACE_SPH_TABLE_H__
#define __DASHMM_LAPLACE_SPH_TABLE_H__


/// \file include/builtins/laplace_sph_table.h
/// \brief Declaration of precomputed tables for LaplaceSPH


#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <vector>
#include "dashmm/types.h"

namespace dashmm {


struct laplace_cmp {
  bool operator()(const double &a, const double &b) const {
    // The smallest gap between key values in the Laplace rotation matrix map is
    // 0.01. The operator compares the first 6 digits and that should be
    // enough.
    double aa = floor(a * 1000000.0) / 100000.0;
    double bb = floor(b * 1000000.0) / 100000.0;
    if (aa == bb)
      return false;
    return aa < bb;
  }
};


using laplace_map_t = std::map<double, double *, laplace_cmp>;

class LaplaceSPHTable {
 public:
  LaplaceSPHTable(int n_digits);
  ~LaplaceSPHTable();

  int p() const {return p_;}
  int s() const {return s_;}
  int nexp() const {return nexp_;}
  const double *sqf() const {return sqf_;}
  const double *sqbinom() const {return sqbinom_;}
  const double *dmat_plus(double v) const {return dmat_plus_->at(v);}
  const double *dmat_minus(double v) const {return dmat_minus_->at(v);}
  const double *lambda() const {return lambda_;}
  const double *weight() const {return weight_;}
  const dcomplex_t *xs() const {return xs_;}
  const dcomplex_t *ys() const {return ys_;}
  const double *zs() const {return zs_;}
  const double *lambdaknm() const {return lambdaknm_;}
  const dcomplex_t *ealphaj() const {return ealphaj_;}
  const int *m() const {return m_;}
  const int *sm() const {return sm_;}
  const int *f() const {return f_;}
  const int *smf() const {return smf_;}

 private:
  int p_;
  double *sqf_;
  double *sqbinom_;
  laplace_map_t *dmat_plus_;
  laplace_map_t *dmat_minus_;

  int s_; 
  int nexp_; 
  double *lambda_; 
  double *weight_; 
  int *m_; 
  int *sm_; 
  int *f_; 
  int *smf_; 
  dcomplex_t *xs_; 
  dcomplex_t *ys_; 
  double *zs_; 
  double *lambdaknm_; 
  dcomplex_t *ealphaj_; 

  double *generate_sqf();
  double *generate_sqbinom();
  void generate_wigner_dmatrix(laplace_map_t *&dp, laplace_map_t *&dm);
  void generate_dmatrix_of_beta(double beta, double *dp, double *dm);
  dcomplex_t *generate_xs(); 
  dcomplex_t *generate_ys(); 
  double *generate_zs(); 
  double *generate_lambdaknm(); 
  dcomplex_t *generate_ealphaj(); 
};

using uLaplaceSPHTable = std::unique_ptr<LaplaceSPHTable>;
using LaplaceSPHTableIterator = std::map<int, uLaplaceSPHTable>::iterator;


extern std::map<int, uLaplaceSPHTable> builtin_laplace_table_;


void legendre_Plm(int n, double x, double *P);


inline int midx(const int n, const int m) {
  return n * (n + 1) / 2 + m;
}


inline int didx(const int n, const int mp, const int m) {
  return n * (n + 1) * (4 * n - 1) / 6 + mp * (2 * n + 1) + n + m;
}


inline double pow_m1(const int m) {
  return (m % 2 ? -1.0 : 1.0);
}


LaplaceSPHTableIterator get_or_add_laplace_sph_table(int n_digits);

enum MergedList {
  uall = 0, ///< +z direction list for all boxes
  u1234 = 1, ///< +z direction list for boxes 1, 2, 3, 4
  nall = 2, ///< +y direction list for all boxes
  n1256 = 3, ///< +y direction list for boxes 1, 2, 5, 6
  n12 = 4, ///< +y direction list for boxes 1, 2 
  n56 = 5, ///< +y direction list for boxes 5, 6
  eall = 6, ///< +x direction list for all boxes 
  e1357 = 7, ///< +x direction list for boxes 1, 3, 5, 7
  e13 = 8, ///< +x direction list for boxes 1, 3
  e57 = 9, ///< +x direction list for boxes 5, 7
  e1 = 10, ///< +x direction list for box 1
  e3 = 11, ///< +x direction list for box 3
  e5 = 12, ///< +x direction list for box 5
  e7 = 13, ///< +x direction list for box 7
  dall = 14, ///< -z direction list for all boxes
  d5678 = 15, ///< -z direction list for boxes 5, 6, 7, 8 
  sall = 16, ///< -y direction list for all boxes
  s3478 = 17, ///< -y direction list for boxes 3, 4, 7, 8 
  s34 = 18, ///< -y direction list for boxes 3, 4 
  s78 = 19, ///< -y direction list for boxes 7, 8 
  wall = 20, ///< -x direction list for all boxes 
  w2468 = 21, ///< -x direction list for boxes 2, 4, 6, 8
  w24 = 22, ///< -x direction list for boxes 2, 4
  w68 = 23, ///< -x direction list for boxes 6, 8
  w2 = 24, ///< -x direction list for box 2
  w4 = 25, ///< -x direction list for box 4
  w6 = 26, ///< -x direction list for box 6
  w8 = 27, ///< -x direction list for box 8
}; 


} // namespace dashmm


#endif // __DASHMM_LAPLACE_SPH_TABLE_H__
