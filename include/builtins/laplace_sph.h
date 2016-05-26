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


#ifndef __DASHMM_LAPLACE_SPH_EXPANSION_H__
#define __DASHMM_LAPLACE_SPH_EXPANSION_H__


/// \file include/builtins/laplace_sph.h
/// \brief Declaration of LaplaceSPH


#include <cassert>
#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <vector>

#include "dashmm/index.h"
#include "builtins/laplace_sph_table.h"
#include "dashmm/point.h"
#include "dashmm/types.h"


namespace dashmm {


struct LaplaceSPHData {
  int reserved;
  int n_digits;
  Point center;
  dcomplex_t expansion[];
};


/// Laplace kernel Spherical Harmonic expansion
///
/// This expansion is of the Laplace Kernel about the center of the node
/// containing the represented sources. The kernel does not include any
/// scaling for physical constants, and so the user will need to multiply
/// results of this expansion by the relevant factors (including a minus sign
/// if needed).
///
/// This expansion is most relevant for interactions that have both signs
/// of the charge.
///
/// This class is a template with parameters for the source and target
/// types.
///
/// Source must define a double valued 'charge' member to be used with
/// LaplaceCOM. Target must define a std::complex<double> valued 'phi' member
/// to be used with LaplaceCOM.
template <typename Source, typename Target>
class LaplaceSPH {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = LaplaceSPH<Source, Target>;

  LaplaceSPH(Point center, int n_digits) {
    LaplaceSPHTableIterator entry = builtin_laplace_table_.find(n_digits);
    assert(entry != builtin_laplace_table_.end());
    uLaplaceSPHTable &table = entry->second;
    int p = table->p();
    int n_terms = (p + 1) * (p + 2) / 2;
    bytes_ = sizeof(LaplaceSPHData) + sizeof(dcomplex_t) * n_terms;
    data_ = reinterpret_cast<LaplaceSPHData *>(new char [bytes_]);
    assert(valid());
    data_->n_digits = n_digits;
    data_->center = center;
    n_digits_ = n_digits;
    for (int i = 0; i < n_terms; ++i)
      data_->expansion[i] = 0;
  }

  LaplaceSPH(void *ptr, size_t bytes, int n_digits)
      : n_digits_{n_digits} {
    data_ = static_cast<LaplaceSPHData *>(ptr);
    bytes_ = bytes;
    if (data_)
      data_->n_digits = n_digits;
  }

  ~LaplaceSPH() {
    if (valid()) {
      delete [] data_;
      data_ = nullptr;
    }
  }

  void *release() {
    LaplaceSPHData *retval = data_;
    data_ = nullptr;
    return retval;
  }

  size_t bytes() const {return bytes_;}

  bool valid() const {return data_ != nullptr;}

  int accuracy() const {return n_digits_;}

  size_t size() const {
    uLaplaceSPHTable &table = builtin_laplace_table_.at(data_->n_digits);
    int p = table->p();
    return (p + 1) * (p + 2) / 2;
  }

  Point center() const {
    assert(valid());
    return data_->center;
  }

  dcomplex_t term(size_t i) const {
    return data_->expansion[i];
  }

  // TODO: This needs to be updated to return an expansion
  std::unique_ptr<expansion_t> S_to_M(Point center, Source *first, Source *last,
                                      double scale) const {
    data_->center = center;
    dcomplex_t *expansion = &data_->expansion[0];
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
    int p = table->p();
    const double *sqf = table->sqf();

    double *legendre = new double[(p + 1) * (p + 2) / 2];
    double *powers_r = new double[p + 1];
    dcomplex_t *powers_ephi = new dcomplex_t[p + 1];
    powers_r[0] = 1.0;
    powers_ephi[0] = dcomplex_t{1.0, 0.0};

    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, center);
      double q = i->charge;
      double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
      double r = dist.norm();

      double ctheta = (r <= 1e-14 ? 1.0 : dist.z() / r);

      // Compute exp(-i * phi) for the azimuthal angle phi
      dcomplex_t ephi = (proj / r <= 1e-14 ? dcomplex_t{1.0, 0.0} :
                         dcomplex_t{dist.x() / proj, -dist.y() / proj});

      // Compute powers of r
      r *= scale;
      for (int j = 1; j <= p; ++j) {
        powers_r[j] = powers_r[j - 1] * r;
      }

      // Compute powers of exp(-i * phi)
      for (int j = 1; j <= p; ++j) {
        powers_ephi[j] = powers_ephi[j - 1] * ephi;
      }

      // Compute multipole expansion M_n^m
      legendre_Plm(p, ctheta, legendre);
      for (int m = 0; m <= p; ++m) {
        for (int n = m; n <= p; ++n) {
          expansion[midx(n, m)] += q * powers_r[n] * powers_ephi[m] *
                                   legendre[midx(n, m)] *
                                   sqf[n - m] / sqf[n + m];
        }
      }
    }

    delete [] legendre;
    delete [] powers_r;
    delete [] powers_ephi;

    // TODO: fix this
    return std::unique_ptr<expansion_t>{nullptr};
  }

  std::unique_ptr<expansion_t> S_to_L(Point center, Source *first,
                                      Source *last, double scale) const {
    expansion_t *retval{new expansion_t{center, n_digits_}};
    dcomplex_t *expansion = &retval->data_->expansion[0];
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
    int p = table->p();
    const double *sqf = table->sqf();

    double *legendre = new double[(p + 1) * (p + 2) / 2];
    double *powers_r = new double[p + 1];
    dcomplex_t *powers_ephi = new dcomplex_t[p + 1];
    powers_ephi[0] = dcomplex_t{1.0, 0.0};

    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, center);
      double q = i->charge;
      double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
      double r = dist.norm();

      // Compute cosine of the polar angle theta
      double ctheta = (r <= 1e-14 ? 1.0 : dist.z() / r);

      // Compute exp(-i * phi) for the azimuthal angle phi
      dcomplex_t ephi = (proj / r <= 1e-14 ? dcomplex_t{1.0, 0.0} :
                         dcomplex_t{dist.x() / proj, -dist.y() / proj});

      // Compute powers of 1 / r
      powers_r[0] = 1.0 / r;
      r *= scale;
      for (int j = 1; j <= p; ++j) {
        powers_r[j] = powers_r[j - 1] / r;
      }

      // Compute powers of exp(-i * phi)
      for (int j = 1; j <= p; ++j) {
        powers_ephi[j] = powers_ephi[j - 1] * ephi;
      }

      // compute local expansion L_n^m
      legendre_Plm(p, ctheta, legendre);
      for (int m = 0; m <= p; ++m) {
        for (int n = m; n <= p; ++n) {
          expansion[midx(n, m)] += q * powers_r[n] * powers_ephi[m] *
                                   legendre[midx(n, m)] *
                                   sqf[n - m] / sqf[n + m];
        }
      }
    }

    delete [] legendre;
    delete [] powers_r;
    delete [] powers_ephi;

    return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> M_to_M(int from_child,
                                      double s_size) const {
    // The function is called on th expansion of the child box and
    // s_size is the child box's size.
    double h = s_size / 2;
    double px = data_->center.x() + (from_child % 2 == 0 ? h : -h);
    double py = data_->center.y() + (from_child % 4 <= 1 ? h : -h);
    double pz = data_->center.z() + (from_child < 4 ? h : -h);

    expansion_t *retval{new expansion_t{Point{px, py, pz}, n_digits_}};

    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
    int p = table->p();
    const double *sqbinom = table->sqbinom();

    // Get precomputed Wigner d-matrix for rotation about the y-axis
    const double *d1 = (from_child < 4 ?
                        table->dmat_plus(1.0 / sqrt(3.0)) :
                        table->dmat_plus(-1.0 / sqrt(3.0)));
    const double *d2 = (from_child < 4 ?
                        table->dmat_minus(1.0 / sqrt(3.0)) :
                        table->dmat_minus(-1.0 / sqrt(3.0)));

    // Shift distance along the z-axis, combined with Y_n^0(pi, 0)
    const double rho = -sqrt(3) / 2;

    // Compute powers of rho
    double *powers_rho = new double[p + 1];
    powers_rho[0] = 1.0;
    for (int i = 1; i <= p; ++i) {
      powers_rho[i] = powers_rho[i - 1] * rho;
    }

    dcomplex_t *W1 = &retval->data_->expansion[0];
    dcomplex_t *W2 = new dcomplex_t[(p + 1) * (p + 2) / 2];

    // Table of rotation angle about the z-axis, as an integer multiple of pi/4
    const int tab_alpha[8] = {1, 3, 7, 5, 1, 3, 7, 5};
    // Get rotation angle about the z-axis
    double alpha = tab_alpha[from_child] * M_PI_4;

    // Rotate the multipole expansion of the child box about z-axis
    dcomplex_t *M = &data_->expansion[0];
    retval->rotate_sph_z(M, alpha, W1);

    // Rotate the previous result further about the y-axis
    retval->rotate_sph_y(W1, d1, W2);

    // Offset to operate multipole expansion
    int offset = 0;

    // Shift along the z-axis by a distance of rho, write result in W1
    for (int n = 0; n <= p; ++n) {
      for (int m = 0; m <= n; ++m) {
        W1[offset] = W2[offset];
        for (int k = 1; k <= n - m; ++k) {
          W1[offset] += W2[midx(n - k, m)] * powers_rho[k] *
                        sqbinom[midx(n - m, k)] * sqbinom[midx(n + m, k)];
        }
        offset++;
      }
    }

    // Reverse rotate the shifted harmonic expansion about the y-axis
    retval->rotate_sph_y(W1, d2, W2);

    // Reverse rotate the previous result further about the z-axis
    retval->rotate_sph_z(W2, -alpha, W1);

    double scale = 1;
    offset = 0;
    for (int n = 0; n <= p; ++n) {
      for (int m = 0; m <= n; ++m) {
        W1[offset++] *= scale;
      }
      scale /= 2;
    }

    delete [] W2;
    delete [] powers_rho;
    return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> M_to_L(Index s_index, double s_size,
                                      Index t_index) const {
    int t2s_x = s_index.x() - t_index.x();
    int t2s_y = s_index.y() - t_index.y();
    int t2s_z = s_index.z() - t_index.z();
    double tx = data_->center.x() - t2s_x * s_size;
    double ty = data_->center.y() - t2s_y * s_size;
    double tz = data_->center.z() - t2s_z * s_size;

    expansion_t *retval{new expansion_t{Point{tx, ty, tz}, n_digits_}};
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
    int p = table->p();


    // Shifting distance
    double rho = sqrt(t2s_x * t2s_x + t2s_y * t2s_y + t2s_z * t2s_z);

    // Compute powers of rho
    double *powers_rho = new double[p * 2 + 1];
    powers_rho[0] = 1.0 / rho;
    for (int i = 1; i <= p * 2; i++) {
      powers_rho[i] = powers_rho[i - 1] / rho;
    }

    // Scaling factor
    double scale = 1.0 / s_size;

    // Temporary space to hold rotated spherical harmonic
    dcomplex_t *W1 = &retval->data_->expansion[0];
    dcomplex_t *W2 = new dcomplex_t[(p + 1) * (p + 2) / 2];

    // Compute the projection of t2s on the x-y plane
    const double proj = sqrt(t2s_x * t2s_x + t2s_y * t2s_y);

    // Handle of the multipole expansion
    dcomplex_t *M = &data_->expansion[0];

    if (proj < 1e-14) {
      if (t2s_z > 0) {
        retval->M_to_L_zp(M, powers_rho, scale, W1);
      } else {
        retval->M_to_L_zm(M, powers_rho, scale, W1);
      }
    } else {
      // azimuthal angle
      double beta = acos(t2s_x / proj);
      if (t2s_y < 0) {
        beta = 2 * M_PI - beta;
      }

      // Get precomputed Wigner d-matrix for rotation about y-axis
      const double *d1 = table->dmat_plus(t2s_z / rho);
      const double *d2 = table->dmat_minus(t2s_z / rho);
      retval->rotate_sph_z(M, beta, W1);
      retval->rotate_sph_y(W1, d1, W2);
      retval->M_to_L_zp(W2, powers_rho, scale, W1);
      retval->rotate_sph_y(W1, d2, W2);
      retval->rotate_sph_z(W2, -beta, W1);
    }

    delete [] W2;
    delete [] powers_rho;
    return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> L_to_L(int to_child, double t_size) const {
    // The function is called on the parent box and t_size is its child size
    double h = t_size / 2;
    double cx = data_->center.x() + (to_child % 2 == 0 ? -h : h);
    double cy = data_->center.y() + (to_child % 4 <= 1 ? -h : h);
    double cz = data_->center.z() + (to_child < 4 ? -h : h);

    expansion_t *retval{new expansion_t{Point{cx, cy, cz}, n_digits_}};

    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
    int p = table->p();
    const double *sqbinom = table->sqbinom();

    // Get precomputed Wigner d-matrix for rotation about the y-axis
    const double *d1 = (to_child < 4 ?
                        table->dmat_plus(1.0 / sqrt(3.0)) :
                        table->dmat_plus(-1.0 / sqrt(3.0)));
    const double *d2 = (to_child < 4 ?
                        table->dmat_minus(1.0 / sqrt(3.0)) :
                        table->dmat_minus(-1.0 / sqrt(3.0)));

    // Shift distance along the z-axis, combined with Y_n^0(pi, 0)
    const double rho = -sqrt(3) / 4;

    // Compute powers of rho
    double *powers_rho = new double[p + 1];
    powers_rho[0] = 1.0;
    for (int i = 1; i <= p; ++i) {
      powers_rho[i] = powers_rho[i - 1] * rho;
    }

    dcomplex_t *W1 = &retval->data_->expansion[0];
    dcomplex_t *W2 = new dcomplex_t[(p + 1) * (p + 2) / 2];

    // Table of rotation angle about the z-axis as an integer multiple of pi / 4
    const int tab_alpha[8] = {1, 3, 7, 5, 1, 3, 7, 5};
    // Get rotation angle about the z-axis
    double alpha = tab_alpha[to_child] * M_PI_4;

    // Rotate the local expansion of the parent box about z-axis
    dcomplex_t *L = &data_->expansion[0];
    retval->rotate_sph_z(L, alpha, W1);

    // Rotate the previous result further about the y-axis
    retval->rotate_sph_y(W1, d1, W2);

    // Offset to operate local expansion
    int offset = 0;

    // Shift along the z-axis by a distance of rho, write result in W1
    for (int n = 0; n <= p; ++n) {
      for (int m = 0; m <= n; ++m) {
        W1[offset] = W2[offset];
        for (int k = 1; k <= p - n; k++) {
          W1[offset] += W2[midx(n + k, m)] * powers_rho[k] *
                        sqbinom[midx(n + k - m, k)] *
                        sqbinom[midx(n + k + m, k)];
        }
        offset++;
      }
    }

    // Reverse rotate the shifted harmonic expansion about the y-axis
    retval->rotate_sph_y(W1, d2, W2);

    // Reverse rotate the previous result further about the z-axis
    retval->rotate_sph_z(W2, -alpha, W1);

    double scale = 1;
    offset = 0;
    for (int n = 0; n <= p; ++n) {
      for (int m = 0; m <= n; ++m) {
        W1[offset++] *= scale;
      }
      scale /= 2;
    }

    delete [] W2;
    delete [] powers_rho;
    return std::unique_ptr<expansion_t>{retval};
  }

  void M_to_T(Target *first, Target *last, double scale) const {
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
    int p = table->p();
    const double *sqf = table->sqf();

    double *legendre = new double[(p + 1) * (p + 2) / 2];
    double *powers_r = new double[p + 1];
    dcomplex_t *powers_ephi = new dcomplex_t[p + 1];
    powers_ephi[0] = dcomplex_t{1.0, 0.0};
    dcomplex_t *M = &data_->expansion[0];

    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, data_->center);
      dcomplex_t potential{0.0, 0.0};
      double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
      double r = dist.norm();

      // Compute cosine of the polar angle theta
      double ctheta = (r <= 1e-14 ? 1.0 : dist.z() / r);

      // Compute exp(i * phi) for the azimuthal angle phi
      dcomplex_t ephi = (proj / r <= 1e-14 ? dcomplex_t{1.0, 0.0} :
                         dcomplex_t{dist.x() / proj, dist.y() / proj});

      // Compute powers of 1 / r
      powers_r[0] = 1.0 / r;
      r *= scale;
      for (int j = 1; j <= p; ++j) {
        powers_r[j] = powers_r[j - 1] / r;
      }

      // Compute powers of exp(i * phi)
      for (int j = 1; j <= p; ++j) {
        powers_ephi[j] = powers_ephi[j - 1] * ephi;
      }

      // Evaluate the multipole expansion M_n^0
      legendre_Plm(p, ctheta, legendre);
      for (int n = 0; n <= p; ++n) {
        potential += M[midx(n, 0)] * powers_r[n] * legendre[midx(n, 0)];
      }

      // Evaluate the multipole expansions M_n^m, where m = 1, ..., p
      for (int m = 1; m <= p; ++m) {
        for (int n = m; n <= p; ++n) {
          potential += 2.0 * real(M[midx(n, m)] * powers_ephi[m]) *
                       powers_r[n] * legendre[midx(n, m)] *
                       sqf[n - m] / sqf[n + m];
        }
      }

      i->phi += potential;
    }

    delete [] powers_r;
    delete [] powers_ephi;
    delete [] legendre;
  }

  void L_to_T(Target *first, Target *last, double scale) const {
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
    int p = table->p();
    const double *sqf = table->sqf();

    double *legendre = new double[(p + 1) * (p + 2) / 2];
    double *powers_r = new double[p + 1];
    dcomplex_t *powers_ephi = new dcomplex_t[p + 1];
    powers_r[0] = 1.0;
    powers_ephi[0] = dcomplex_t{1.0, 0.0};

    dcomplex_t *L = &data_->expansion[0];

    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, data_->center);
      dcomplex_t potential{0.0, 0.0};
      double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
      double r = dist.norm();

      // Compute cosine of the polar angle theta
      double ctheta = (r <= 1e-14 ? 1.0 : dist.z() / r);

      // Compute exp(i * phi) for the azimuthal angle phi
      dcomplex_t ephi = (proj / r <= 1e-14 ? dcomplex_t{1.0, 0.0} :
                         dcomplex_t{dist.x() / proj, dist.y() / proj});

      // Compute powers of r
      r *= scale;
      for (int j = 1; j <= p; ++j) {
        powers_r[j] = powers_r[j - 1] * r;
      }

      // Compute powers of exp(i * phi)
      for (int j = 1; j <= p; ++j) {
        powers_ephi[j] = powers_ephi[j - 1] * ephi;
      }

      // Evaluate the local expansion L_n^0
      legendre_Plm(p, ctheta, legendre);
      for (int n = 0; n <= p; ++n) {
        potential += L[midx(n, 0)] * powers_r[n] * legendre[midx(n, 0)];
      }

      // Evaluate the local expansions L_n^m, where m = 1, ..., p
      for (int m = 1; m <= p; ++m) {
        for (int n = m; n <= p; ++n) {
          potential += 2.0 * real(L[midx(n, m)] * powers_ephi[m]) *
                       powers_r[n] * legendre[midx(n, m)] *
                       sqf[n - m] / sqf[n + m];
        }
      }

      i->phi += potential;
    }

    delete [] powers_r;
    delete [] powers_ephi;
    delete [] legendre;
  }

  void S_to_T(Source *s_first, Source *s_last,
              Target *t_first, Target *t_last) const {
    for (auto i = t_first; i != t_last; ++i) {
      dcomplex_t potential{0.0, 0.0};
      for (auto j = s_first; j != s_last; ++j) {
        Point s2t = point_sub(i->position, j->position);
        double dist = s2t.norm();
        if (dist > 0)
          potential += j->charge / dist;
      }
      i->phi += potential;
    }
  }

  void add_expansion(const expansion_t *temp1) {
    dcomplex_t *expansion = &data_->expansion[0];
    for (size_t i = 0; i < temp1->size(); ++i)
      expansion[i] += temp1->term(i);
  }

 private:
  LaplaceSPHData *data_;
  size_t bytes_;
  int n_digits_;

  void rotate_sph_z(const dcomplex_t *M, double alpha, dcomplex_t *MR) {
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
    int p = table->p();

    // Compute exp(i * alpha)
    dcomplex_t ealpha{cos(alpha), sin(alpha)};

    // Compute powers of exp(i * alpha)
    dcomplex_t *powers_ealpha = new dcomplex_t[p + 1];
    powers_ealpha[0] = dcomplex_t{1.0, 0.0};
    for (int j = 1; j <= p; ++j) {
      powers_ealpha[j] = powers_ealpha[j - 1] * ealpha;
    }

    int offset = 0;
    for (int n = 0; n <= p; n++) {
      for (int m = 0; m <= n; m++) {
        MR[offset] = M[offset] * powers_ealpha[m];
        offset++;
      }
    }

    delete [] powers_ealpha;
  }

  void rotate_sph_y(const dcomplex_t *M, const double *d, dcomplex_t *MR) {
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
    int p = table->p();

    int offset = 0;
    for (int n = 0; n <= p; ++n) {
      int power_mp = 1;
      for (int mp = 0; mp <= n; ++mp) {
        // Retrieve address of wigner d-matrix entry d_n^{mp, 0}
        const double *coeff = &d[didx(n, mp, 0)];
        // Get address of original harmonic expansion M_n^0
        const dcomplex_t *Mn = &M[midx(n, 0)];
        // Compute rotated spherical harmonic M_n^mp
        MR[offset] = Mn[0] * coeff[0];
        double power_m = -1;
        for (int m = 1; m <= n; ++m) {
          MR[offset] += (Mn[m] * power_m * coeff[m] + conj(Mn[m]) * coeff[-m]);
          power_m = -power_m;
        }
        MR[offset++] *= power_mp;
        power_mp = -power_mp;
      }
    }
  }

  void M_to_L_zp(const dcomplex_t *M, const double *rho, double scale,
                 dcomplex_t *L) {
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
    int p = table->p();
    const double *sqbinom = table->sqbinom();

    int offset = 0;
    for (int j = 0; j <= p; ++j) {
      for (int k = 0; k <= j; ++k) {
        L[offset] = 0;
        for (int n = k; n <= p; ++n) {
          L[offset] += M[midx(n, k)] * pow_m1(n + k) * rho[j + n] *
                       sqbinom[midx(n + j, n - k)] * sqbinom[midx(n + j, n + k)];
        }
        L[offset] *= scale;
        offset++;
      }
    }
  }

  void M_to_L_zm(const dcomplex_t *M, const double *rho, double scale,
                 dcomplex_t *L) {
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits_);
    int p = table->p();
    const double *sqbinom = table->sqbinom();

    int offset = 0;
    for (int j = 0; j <= p; ++j) {
      for (int k = 0; k <= j; ++k) {
        L[offset] = 0;
        for (int n = k; n <= p; ++n) {
          L[offset] += M[midx(n, k)] * pow_m1(k + j) * rho[j + n] *
                       sqbinom[midx(n + j, n - k)] * sqbinom[midx(n + j, n + k)];
        }
        L[offset] *= scale;
        offset++;
      }
    }
  }
};


/// Precompute some required data for using the LaplaceSPH expansion
///
/// This routine must be called before an object of type LaplaceSPH is
/// constucted with the given accuracy.
///
/// \param n_digits - the accuracy parameter for which to pregenerate values
void laplace_sph_precompute(int n_digits);


} // namespace dashmm

#endif // __DASHMM_LAPLACE_SPH_EXPANSION_H__
