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


#ifndef __DASHMM_LAPLACE_TABLE_H__
#define __DASHMM_LAPLACE_TABLE_H__


/// \file
/// \brief Declaration of precomputed tables for Laplace


#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <vector>
#include "dashmm/types.h"
#include "builtins/special_function.h"
#include "builtins/builtin_table.h"

namespace dashmm {

class LaplaceTable {
 public:
  LaplaceTable(int n_digits, double size);
  ~LaplaceTable();

  int p() const {return p_;}
  int s() const {return s_;}
  int nexp() const {return nexp_;}
  double scale(int lev) const {return scale_ * pow(2, lev);}
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
  double scale_; // scaling factor of level 0 to normalize box size to 1
  double *sqf_;
  double *sqbinom_;
  builtin_map_t *dmat_plus_;
  builtin_map_t *dmat_minus_;

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

  void generate_sqf();
  void generate_sqbinom();
  void generate_wigner_dmatrix();
  void generate_dmatrix_of_beta(double beta, double *dp, double *dm);
  void generate_xs();
  void generate_ys();
  void generate_zs();
  void generate_lambdaknm();
  void generate_ealphaj();
};

extern std::unique_ptr<LaplaceTable> builtin_laplace_table_;

void update_laplace_table(int n_digits, double size);

} // namespace dashmm


#endif // __DASHMM_LAPLACE_TABLE_H__
