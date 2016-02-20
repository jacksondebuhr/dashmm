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


/// \file include/laplace_sph_table.h
/// \brief Declaration of precomputed tables for LaplaceSPH


#include <cmath>
#include <complex>
#include <map>
#include <vector>


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
  const double *sqf() const {return sqf_;}
  const double *sqbinom() const {return sqbinom_;}
  const double *dmat_plus(double v) const {return dmat_plus_->at(v);}
  const double *dmat_minus(double v) const {return dmat_minus_->at(v);}

 private:
  int p_;
  double *sqf_;
  double *sqbinom_;
  laplace_map_t *dmat_plus_;
  laplace_map_t *dmat_minus_;

  double *generate_sqf();
  double *generate_sqbinom();
  void generate_wigner_dmatrix(laplace_map_t *&dp, laplace_map_t *&dm);
  void generate_dmatrix_of_beta(double beta, double *dp, double *dm);
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


} // namespace dashmm


#endif // __DASHMM_LAPLACE_SPH_TABLE_H__
