// =============================================================================
//  This file is part of:
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  DASHMM is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  DASHMM is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with DASHMM. If not, see <http://www.gnu.org/licenses/>.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================


#ifndef __DASHMM_YUKAWA_TABLE_H__
#define __DASHMM_YUKAWA_TABLE_H__


/// \file
/// \brief Declaration of precomputed tables for Yukawa


#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <vector>
#include "dashmm/types.h"
#include "builtins/special_function.h"
#include "builtins/builtin_table.h"


namespace dashmm {


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
  builtin_map_t *dmat_plus_;
  builtin_map_t *dmat_minus_;
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

extern std::unique_ptr<YukawaTable> builtin_yukawa_table_;

void update_yukawa_table(int n_digits, double size, double lambda);


} // namespace dashmm

#endif // __DASHMM_YUKAWA_TABLE_H__
