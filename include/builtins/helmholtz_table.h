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

#ifndef __DASHMM_HELMHOLTZ_TABLE_H__
#define __DASHMM_HELMHOLTZ_TABLE_H__

/// \file
/// \brief Declaration of precomputed table for Helmholtz

#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <vector>
#include "dashmm/types.h"
#include "builtins/special_function.h"
#include "builtins/builtin_table.h"

namespace dashmm {

class HelmholtzTable {
public:
  HelmholtzTable(int n_digits, double size, double omega);
  ~HelmholtzTable();
  static const int maxlev;
  int p() const {return p_;}
  int s_e() const {return s_e_;}
  int s_p() const {return s_p_;}
  double scale(int lev) const {return scale_ / pow(2, lev);}
  double omega() const {return omega_;}
  const double *sqf() const {return sqf_;}
  const double *dmat_plus(double v) const {return dmat_plus_->at(v);}
  const double *dmat_minus(double v) const {return dmat_minus_->at(v);}
  const double *m2m(double scale) const {
    int offset = level(scale) - 3;
    assert(offset >= 0);
    return &m2m_[(p_ + 1) * (p_ + 1) * (p_ + 2) / 2 * offset];
  }
  const double *l2l(double scale) const {
    int offset = level(scale) - 2;
    assert(offset >= 0 && offset <= maxlev - 3);
    return &l2l_[(p_ + 1) * (p_ + 1) * (p_ + 2) / 2 * offset];
  }
  const double *x_e() const {return x_e_;}
  const double *w_e() const {return w_e_;}
  const int *m_e(double scale) const {
    return &m_e_[s_e_ * level(scale)];
  }
  const int *sm_e(double scale) const {
    return &sm_e_[(s_e_ + 1) * level(scale)];
  }
  int n_e(double scale) const {return n_e_[level(scale)];}
  const int *f_e() const {return f_e_;}
  const int *smf_e(double scale) const {
    return &smf_e_[(s_e_ + 1) * level(scale)];
  }
  const dcomplex_t *ealphaj_e(double scale) const {
    return ealphaj_e_[level(scale)];
  }
  const dcomplex_t *xs_e(double scale) const {
    return xs_e_[level(scale)];
  }
  const dcomplex_t *ys_e(double scale) const {
    return ys_e_[level(scale)];
  }
  const double *zs_e(double scale) const {
    return zs_e_[level(scale)];
  }

  const double *x_p() const {return x_p_;}
  const double *w_p() const {return w_p_;}
  const int *m_p(double scale) const {
    return &m_p_[s_p_ * level(scale)];
  }
  const int *sm_p(double scale) const {
    return &sm_p_[(s_p_ + 1) * level(scale)];
  }
  int n_p(double scale) const {return n_p_[level(scale)];}
  const int *f_p() const {return f_p_;}
  const int *smf_p(double scale) const {
    return &smf_p_[(s_p_ + 1) * level(scale)];
  }
  const dcomplex_t *ealphaj_p(double scale) const {
    return ealphaj_p_[level(scale)];
  }
  const dcomplex_t *xs_p(double scale) const {
    return xs_p_[level(scale)];
  }
  const dcomplex_t *ys_p(double scale) const {
    return ys_p_[level(scale)];
  }
  const dcomplex_t *zs_p(double scale) const {
    return zs_p_[level(scale)];
  }
  double size(double scale) const {return size_ * scale / scale_;}

  bool update(int n_digits, double size, double omega) const {
    if (n_digits_ != n_digits || size_ != size || omega_ != omega) 
      return true;
    return false;
  }
    
private:
  int p_;
  int s_e_;
  int s_p_;
  int n_digits_; 
  double omega_;
  double size_;
  double scale_;
  double *sqf_;
  builtin_map_t *dmat_plus_;
  builtin_map_t *dmat_minus_;
  double *m2m_;
  double *l2l_;
  double *x_e_;
  double *w_e_;
  int *m_e_;
  int *sm_e_;
  int *n_e_;
  int *f_e_;
  int *smf_e_;
  dcomplex_t **ealphaj_e_;
  dcomplex_t **xs_e_;
  dcomplex_t **ys_e_;
  double **zs_e_;
  double *x_p_;
  double *w_p_;
  int *m_p_;
  int *sm_p_;
  int *n_p_;
  int *f_p_;
  int *smf_p_;
  dcomplex_t **ealphaj_p_;
  dcomplex_t **xs_p_;
  dcomplex_t **ys_p_;
  dcomplex_t **zs_p_;

  int level(double scale) const {return log2(scale_ / scale);}
  void gaussq(int N);
  int imtql2(int N, double *d, double *e, double *z);
  void generate_sqf();
  void generate_scaled_wigner_dmat();
  void generate_scaled_dmat_of_beta(double beta, double *dp, double *dm);
  void generate_m2m();
  void generate_l2l();
};


extern std::unique_ptr<HelmholtzTable> builtin_helmholtz_table_;

void update_helmholtz_table(int n_digits, double size, double omega);

}; // namespace dashmm

#endif
