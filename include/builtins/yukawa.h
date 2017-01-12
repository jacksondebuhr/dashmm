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


#ifndef __DASHMM_YUKAWA_EXPANSION_H__
#define __DASHMM_YUKAWA_EXPANSION_H__


/// \file
/// \brief Declaration of Yukawa


#include <cassert>
#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <vector>

#include "dashmm/index.h"
#include "builtins/yukawa_table.h"
#include "dashmm/point.h"
#include "dashmm/types.h"
#include "dashmm/viewset.h"


namespace dashmm {


/// Yukawa kernel Spherical Harmonic expansion
///
/// This expansion is of the Yukawa Kernel about the center of the node
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
/// Yukawa. Target must define a std::complex<double> valued 'phi' member
/// to be used with Yukawa.
template <typename Source, typename Target>
class Yukawa {
public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Yukawa<Source, Target>;

  Yukawa(Point center, double scale, ExpansionRole role)
    : views_{ViewSet{role, center, scale}} {

    // View size for each spherical harmonic expansion
    int p = builtin_yukawa_table_->p();
    int nsh = (p + 1) * (p + 2) / 2;

    if (role == kSourcePrimary || role == kTargetPrimary) {
      size_t bytes = sizeof(dcomplex_t) * nsh;
      char *data = new char[bytes]();
      views_.add_view(0, bytes, data);
    } else {
      // View size for each exponential expansion at the current scale level
      int nexp = builtin_yukawa_table_->nexp(scale);

      if (role == kSourceIntermediate) {
        size_t bytes = sizeof(dcomplex_t) * nexp;
        for (int i = 0; i < 6; ++i) {
          char *data = new char[bytes]();
          views_.add_view(i, bytes, data);
        }
      } else { // role == kTargetIntermediate
        size_t bytes = sizeof(dcomplex_t) * nexp;
        for (int i = 0; i < 28; ++i) {
          char *data = new char[bytes]();
          views_.add_view(i, bytes, data);
        }
      }
    }
  }

  Yukawa(const ViewSet &views) : views_{views} { }

  ~Yukawa() {
    int count = views_.count();
    if (count) {
      for (int i = 0; i < count; ++i) {
        delete [] views_.view_data(i);
      }
    }
  }

  void release() {views_.clear();}

  bool valid(const ViewSet &view) const {
    // \p view is assumed to be a subset of \p views_ (no range checking
    // performed). The function returns true if and only if each entry in the
    // required subset is associated with some data.
    bool is_valid = true;
    int count = view.count();
    for (int i = 0; i < count; ++i) {
      int idx = view.view_index(i);
      if (views_.view_data(idx) == nullptr) {
        is_valid = false;
        break;
      }
    }
    return is_valid;
  }

  int view_count() const { return views_.count(); }

  // This is likely to be removed from the interface
  void get_views(ViewSet &view) const {}

  ViewSet get_all_views() const {return views_;}

  ExpansionRole role() const {return views_.role();}

  Point center() const {return views_.center();}

  size_t view_size(int view) const {
    return views_.view_bytes(view) / sizeof(dcomplex_t);
  }

  dcomplex_t view_term(int view, size_t i) const {
    dcomplex_t *data = reinterpret_cast<dcomplex_t *>(views_.view_data(view));
    return data[i];
  }

  std::unique_ptr<expansion_t> S_to_M(Point center, Source *first,
                                      Source *last) const {
    double scale = views_.scale();
    expansion_t *retval{new expansion_t{center, scale, kSourcePrimary}};
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    int p = builtin_yukawa_table_->p();
    const double *sqf = builtin_yukawa_table_->sqf();
    double lambda = builtin_yukawa_table_->lambda();

    double *legendre = new double[(p + 1) * (p + 2) / 2];
    dcomplex_t *powers_ephi = new dcomplex_t[p + 1];
    double *bessel = new double[p + 1];

    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, center);
      double q = i->charge;
      double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
      double r = dist.norm();

      // Compute cosine of the polar angle theta
      double ctheta = (r <= 1.0e-14 ? 1.0 : dist.z() / r);

      // Compute exp(-i * phi) for the azimuthal angle phi
      dcomplex_t ephi = (proj / r <= 1.0e-14 ? dcomplex_t{1.0, 0.0} :
                         dcomplex_t{dist.x() / proj, -dist.y() /proj});

      // Compute powers of exp(-i * phi)
      powers_ephi[0] = 1.0;
      for (int j = 1; j <= p; ++j) {
        powers_ephi[j] = powers_ephi[j - 1] * ephi;
      }

      // Compute scaled modified spherical bessel function
      bessel_in_scaled(p, lambda * r, scale, bessel);

      // Compute legendre polynomial
      legendre_Plm(p, ctheta, legendre);

      // Compute multipole expansion M_n^m
      for (int n = 0; n <= p; ++n) {
        for (int m = 0; m <= n; ++m) {
          int idx = midx(n, m);
          M[idx] += q * bessel[n] * sqf[idx] * legendre[idx] * powers_ephi[m];
        }
      }
    }

    delete [] legendre;
    delete [] powers_ephi;
    delete [] bessel;

    return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> S_to_L(Point center, Source *first,
                                      Source *last) const {
    double scale = views_.scale();
    expansion_t *retval{new expansion_t{center, scale, kTargetPrimary}};
    dcomplex_t *L = reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    int p = builtin_yukawa_table_->p();
    const double *sqf = builtin_yukawa_table_->sqf();
    double lambda = builtin_yukawa_table_->lambda();

    double *legendre = new double[(p + 1) * (p + 2) / 2];
    double *bessel = new double[p + 1];
    dcomplex_t *powers_ephi = new dcomplex_t[p + 1];

    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, center);
      double q = i->charge;
      double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
      double r = dist.norm();

      // Compute cosine of the polar angle theta
      double ctheta = (r <= 1.0e-14 ? 1.0 : dist.z() / r);

      // Compute exp(-i * phi) for the azimuthal angle phi
      dcomplex_t ephi = (proj / r <= 1.0e-14 ? dcomplex_t{1.0, 0.0} :
                         dcomplex_t{dist.x() / proj, -dist.y() / proj});

      // Compute powers of exp(-i * phi)
      powers_ephi[0] = 1.0;
      for (int j = 1; j <= p; j++) {
        powers_ephi[j] = powers_ephi[j - 1] * ephi;
      }

      // Compute scaled modified spherical bessel function
      bessel_kn_scaled(p, lambda * r, scale, bessel);

      // Compute legendre polynomial
      legendre_Plm(p, ctheta, legendre);

      // Compute local expansion L_n^m
      for (int n = 0; n <= p; ++n) {
        for (int m = 0; m <= n; ++m) {
          int idx = midx(n, m);
          L[idx] += q * bessel[n] * sqf[idx] * legendre[idx] * powers_ephi[m];
        }
      }
    }

    delete [] legendre;
    delete [] bessel;
    delete [] powers_ephi;
    return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> M_to_M(int from_child) const {
    double scale = views_.scale();
    expansion_t *retval{new expansion_t{Point{0.0, 0.0, 0.0}, 
          0.0, kSourcePrimary}};
    

    /*
    // The function is called on the expansion of the child box and \p s_size is
    // the child box's size.
    double h = s_size / 2;
    Point center = views_.center();
    double px = center.x() + (from_child % 2 == 0 ? h : -h);
    double py = center.y() + (from_child % 4 <= 1 ? h : -h);
    double pz = center.z() + (from_child < 4 ? h : -h);

    expansion_t *retval{new
        expansion_t{Point{px, py, pz}, scale * 2, kSourcePrimary}};
    */
    int p = builtin_yukawa_table_->p();

    // Get precomputed Wigner d-matrix for rotation about the y-axis
    const double *d1 = (from_child < 4 ?
                        builtin_yukawa_table_->dmat_plus(1.0 / sqrt(3.0)) :
                        builtin_yukawa_table_->dmat_plus(-1.0 / sqrt(3.0)));
    const double *d2 = (from_child < 4 ?
                        builtin_yukawa_table_->dmat_minus(1.0 / sqrt(3.0)) :
                        builtin_yukawa_table_->dmat_minus(-1.0 / sqrt(3.0)));

    // Get precomputed coefficients for shifting along z-axis
    const double *coeff = builtin_yukawa_table_->m2m(scale);

    // Table of rotation angle about z-axis, as an integer of multiple of pi / 4
    const int tab_alpha[8] = {1, 3, 7, 5, 1, 3, 7, 5};

    // Get rotation angle
    double alpha = tab_alpha[from_child] * M_PI_4;

    // Get multipole expansion of the child box
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(views_.view_data(0));

    // Temporary space for rotating multipole expansion
    dcomplex_t *W1 =
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    dcomplex_t *W2 = new dcomplex_t[(p + 1) * (p + 2) / 2];

    rotate_sph_z(M, alpha, W1);
    rotate_sph_y(W1, d1, W2);

    for (int n = 0; n <= p; ++n) {
      for (int m = 0; m <= n; ++m) {
        dcomplex_t temp{0.0, 0.0};
        for (int np = m; np <= p; ++np) {
          temp += W2[midx(np, m)] * coeff[sidx(n, m, np, p)];
        }
        W1[midx(n, m)] = temp;
      }
    }

    rotate_sph_y(W1, d2, W2);
    rotate_sph_z(W2, -alpha, W1);

    delete [] W2;
    return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> M_to_L(Index s_index, Index t_index) const {
    return std::unique_ptr<expansion_t>{nullptr};
  }

  std::unique_ptr<expansion_t> L_to_L(int to_child) const {
    double scale = views_.scale();
    expansion_t *retval{new expansion_t{Point{0.0, 0.0, 0.0}, 
          0.0, kTargetPrimary}}; 

    // Table of rotation angle about z-axis, as an integer multiple of pi / 4
    const int tab_alpha[8] = {1, 3, 7, 5, 1, 3, 7, 5};

    // Get rotation angle
    double alpha = tab_alpha[to_child] * M_PI_4;

    int p = builtin_yukawa_table_->p();

    // Get precomputed Wigner d-matrix for rotation about the y-axis
    const double *d1 = (to_child < 4 ?
                        builtin_yukawa_table_->dmat_plus(1.0 / sqrt(3)) :
                        builtin_yukawa_table_->dmat_plus(-1.0 / sqrt(3)));
    const double *d2 = (to_child < 4 ?
                        builtin_yukawa_table_->dmat_minus(1.0 / sqrt(3)) :
                        builtin_yukawa_table_->dmat_minus(-1.0 / sqrt(3)));

    // Get precomputed coefficients for shifting along z-axis
    const double *coeff = builtin_yukawa_table_->l2l(scale);

    // Get local expansion of the parent box
    dcomplex_t *L = reinterpret_cast<dcomplex_t *>(views_.view_data(0));

    // Temporary space for rotating local expansion
    dcomplex_t *W1 =
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    dcomplex_t *W2 = new dcomplex_t[(p + 1) * (p + 2) / 2];

    rotate_sph_z(L, alpha, W1);
    rotate_sph_y(W1, d1, W2);

    for (int n = 0; n <= p; ++n) {
      for (int m = 0; m <= n; ++m) {
        dcomplex_t temp{0.0, 0.0};
        for (int np = m; np <= p; ++np) {
          temp += W2[midx(np, m)] * coeff[sidx(n, m, np, p)];
        }
        W1[midx(n, m)] = temp;
      }
    }

    rotate_sph_y(W1, d2, W2);
    rotate_sph_z(W2, -alpha, W1);

    delete [] W2;
    return std::unique_ptr<expansion_t>{retval};
  }

  void M_to_T(Target *first, Target *last) const {
    int p = builtin_yukawa_table_->p();
    double scale = views_.scale();
    double lambda = builtin_yukawa_table_->lambda();
    double *legendre = new double[(p + 1) * (p + 2) / 2];
    double *bessel = new double[p + 1];
    dcomplex_t *powers_ephi = new dcomplex_t[p + 1];
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(views_.view_data(0));

    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, views_.center());
      dcomplex_t potential{0.0, 0.0};
      double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
      double r = dist.norm();

      // Compute cosine of the polar angle theta
      double ctheta = (r <= 1.0e-14 ? 1.0 : dist.z() / r);

      // Compute exp(i * phi) for the azimuthal angle phi
      dcomplex_t ephi = (proj / r <= 1.0e-14 ? dcomplex_t{1.0, 0.0} :
                         dcomplex_t{dist.x() / proj, dist.y() / proj});

      // Compute powers of exp(i * phi)
      powers_ephi[0] = 1.0;
      for (int j = 1; j <= p; j++) {
        powers_ephi[j] = powers_ephi[j - 1] * ephi;
      }

      // Compute scaled modified spherical bessel function
      bessel_kn_scaled(p, lambda * r, scale, bessel);

      // Compute legendre polynomial
      legendre_Plm(p, ctheta, legendre);

      // Evaluate M_n^0
      for (int n = 0; n <= p; ++n) {
        potential += M[midx(n, 0)] * legendre[midx(n, 0)] * bessel[n];
      }

      // Evaluate M_n^m
      for (int n = 1; n <= p; ++n) {
        for (int m = 1; m <= n; ++m) {
          potential += 2.0 * real(M[midx(n, m)] * powers_ephi[m]) *
            bessel[n] * legendre[midx(n, m)];
        }
      }

      i->phi += potential;
    }

    delete [] legendre;
    delete [] powers_ephi;
    delete [] bessel;
  }

  void L_to_T(Target *first, Target *last) const {
    int p = builtin_yukawa_table_->p();
    double scale = views_.scale();
    double lambda = builtin_yukawa_table_->lambda();
    double *legendre = new double[(p + 1) * (p + 2) / 2];
    double *bessel = new double[p + 1];
    dcomplex_t *powers_ephi = new dcomplex_t[p + 1];
    dcomplex_t *L = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    //double scale = views_.scale();

    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, views_.center());
      dcomplex_t potential{0.0, 0.0};
      double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
      double r = dist.norm();

      // Compute cosine of the polar angle theta
      double ctheta = (r <= 1.0e-14 ? 1.0 : dist.z() / r);

      // Compute exp(i * phi) for the azimuthal angle phi
      dcomplex_t ephi = (proj / r <= 1.0e-14 ? dcomplex_t{1.0, 0.0} :
                         dcomplex_t{dist.x() / proj, dist.y() / proj});

      // Compute powers of exp(i * phi)
      powers_ephi[0] = 1.0;
      for (int j = 1; j <= p; ++j) {
        powers_ephi[j] = powers_ephi[j - 1] * ephi;
      }

      // Compute scaled modified spherical bessel function
      bessel_in_scaled(p, lambda * r, scale, bessel);

      // Compute legendre polynomial
      legendre_Plm(p, ctheta, legendre);

      // Evaluate local expansion L_n^0
      for (int n = 0; n <= p; ++n) {
        potential += L[midx(n, 0)] * legendre[midx(n, 0)] * bessel[n];
      }

      // Evaluate L_n^m
      for (int n = 1; n <= p; ++n) {
        for (int m = 1; m <= n; ++m) {
          int idx = midx(n, m);
          potential += 2.0 * real(L[idx] * powers_ephi[m]) *
            bessel[n] * legendre[idx];
        }
      }

      i->phi += potential;
    }

    delete [] powers_ephi;
    delete [] bessel;
    delete [] legendre;
  }

  void S_to_T(Source *s_first, Source *s_last,
              Target *t_first, Target *t_last) const {
    double lambda = builtin_yukawa_table_->lambda();
    for (auto i = t_first; i != t_last; ++i) {
      dcomplex_t potential{0.0, 0.0};
      for (auto j = s_first; j != s_last; ++j) {
        Point s2t = point_sub(i->position, j->position);
        double dist = lambda * s2t.norm();
        if (dist > 0) {
          potential += j->charge * exp(-dist) / dist;
        }
      }
      i->phi += potential * M_PI_2;
    }
  }

  std::unique_ptr<expansion_t> M_to_I() const {
    double scale = views_.scale();
    expansion_t *retval{new expansion_t{views_.center(),
          scale, kSourceIntermediate}};
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(views_.view_data(0));

    // Addresses of the views
    dcomplex_t *E_px =
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    dcomplex_t *E_mx =
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(1));
    dcomplex_t *E_py =
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(2));
    dcomplex_t *E_my =
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(3));
    dcomplex_t *E_pz =
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(4));
    dcomplex_t *E_mz =
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(5));

    // Addresses of exponential expansions in the positive axis direction
    dcomplex_t *EP[3] = {E_px, E_py, E_pz};
    // Addresses of exponential expansions in the negative axis direction
    dcomplex_t *EM[3] = {E_mx, E_my, E_mz};

    int p = builtin_yukawa_table_->p();
    double ld = builtin_yukawa_table_->lambda() *
      builtin_yukawa_table_->size(scale);

    // Get precomputed Wigner d-matrix
    const double *d1 = builtin_yukawa_table_->dmat_plus(0.0);
    const double *d2 = builtin_yukawa_table_->dmat_minus(0.0);

    // Get number of Gaussian quadrature points
    int s = builtin_yukawa_table_->s();

    // Get Gaussian quadrature points
    const double *x = builtin_yukawa_table_->x();

    // Get number of Fourier modes
    const int *f = builtin_yukawa_table_->f();

    // Get number of trapezoid points
    const int *m = builtin_yukawa_table_->m(scale);
    const int *smf = builtin_yukawa_table_->smf(scale);

    // Get e^{i * m * alpha_j}
    const dcomplex_t *ealphaj = builtin_yukawa_table_->ealphaj(scale);

    // Allocate temporary space to handle x-/y-direction expansion
    dcomplex_t *W1 = new dcomplex_t[(p + 1) * (p + 2) / 2];
    dcomplex_t *W2 = new dcomplex_t[(p + 1) * (p + 2) / 2];

    // Setup y-direction
    rotate_sph_z(M, -M_PI / 2, W1);
    rotate_sph_y(W1, d2, W2);

    // Setup x-direction
    rotate_sph_y(M, d1, W1);

    // Addresses of the spherical harmonic expansions
    const dcomplex_t *SH[3] = {W1, W2, M};

    double *legendre = new double[(p + 1) * (p + 2) / 2];

    for (int dir = 0; dir <=2; ++dir) {
      int offset = 0;
      for (int k = 0; k < s; ++k) {
        legendre_Plm_gt1_scaled(p, 1 + x[k] / ld, scale, legendre);

        // Handle M_n^m where n is even
        dcomplex_t *z1 = new dcomplex_t[f[k] + 1];
        // Handle M_n^m where n is odd
        dcomplex_t *z2 = new dcomplex_t[f[k] + 1];

        // Process M_n^0 terms
        z1[0] = 0;
        z2[0] = 0;
        for (int n = 0; n <= p; n += 2) {
          z1[0] += SH[dir][midx(n, 0)] * legendre[midx(n, 0)];
        }
        for (int n = 1; n <= p; n += 2) {
          z2[0] += SH[dir][midx(n, 0)] * legendre[midx(n, 0)];
        }

        // Process M_n^m, where m is odd
        for (int m = 1; m <= f[k]; m += 2) {
          z1[m] = 0;
          z2[m] = 0;
          for (int n = m; n <= p; n += 2) {
            z2[m] += SH[dir][midx(n, m)] * legendre[midx(n, m)];
          }
          for (int n = m + 1; n <= p; n += 2) {
            z1[m] += SH[dir][midx(n, m)] * legendre[midx(n, m)];
          }
        }

        // Process M_n^m, where m is even
        for (int m = 2; m <= f[k]; m += 2) {
          z1[m] = 0;
          z2[m] = 0;
          for (int n = m; n <= p; n += 2) {
            z1[m] += SH[dir][midx(n, m)] * legendre[midx(n, m)];
          }
          for (int n = m + 1; n <= p; n += 2) {
            z2[m] += SH[dir][midx(n, m)] * legendre[midx(n, m)];
          }
        }

        // Compute W(k, j)
        for (int j = 1; j <= m[k] / 2; ++j) {
          dcomplex_t up{z1[0] + z2[0]};
          dcomplex_t dn{z1[0] - z2[0]};
          dcomplex_t power_I{0.0, 1.0};
          for (int m = 1; m <= f[k]; ++m) {
            int idx = smf[k] + (j - 1) * f[k] + m - 1;
            up += 2 * real(ealphaj[idx] * (z1[m] + z2[m])) * power_I;
            dn += 2 * real(ealphaj[idx] * (z1[m] - z2[m])) * power_I;
            power_I *= dcomplex_t{0.0, 1.0};
          }
          EP[dir][offset] = up;
          EM[dir][offset] = dn;
          offset++;
        }

        delete [] z1;
        delete [] z2;
      }
    }

    delete [] W1;
    delete [] W2;
    delete [] legendre;
    return std::unique_ptr<expansion_t>(retval);
  }

  std::unique_ptr<expansion_t> I_to_I(Index s_index, double s_size,
                                      Index t_index) const {
    // t_index is the index of the parent node on the target side

    // Compute index offsets between the current source node and the 1st child
    // of the parent node
    int dx = s_index.x() - t_index.x() * 2;
    int dy = s_index.y() - t_index.y() * 2;
    int dz = s_index.z() - t_index.z() * 2;

    // Compute center of the parent node
    Point center = views_.center();
    double px = center.x() + (dx + 0.5) * s_size;
    double py = center.y() + (dy + 0.5) * s_size;
    double pz = center.z() + (dz + 0.5) * s_size;

    // Exponential expansions on the source side
    double scale = views_.scale();
    int nexp = builtin_yukawa_table_->nexp(scale);
    const dcomplex_t *S_px =
      reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    const dcomplex_t *S_mx =
      reinterpret_cast<dcomplex_t *>(views_.view_data(1));
    const dcomplex_t *S_py =
      reinterpret_cast<dcomplex_t *>(views_.view_data(2));
    const dcomplex_t *S_my =
      reinterpret_cast<dcomplex_t *>(views_.view_data(3));
    const dcomplex_t *S_pz =
      reinterpret_cast<dcomplex_t *>(views_.view_data(4));
    const dcomplex_t *S_mz =
      reinterpret_cast<dcomplex_t *>(views_.view_data(5));

    ViewSet views{kTargetIntermediate, Point{px, py, pz}, 2 * scale};

    // Each S is going to generate between 1 and 3 views of the exponential
    // expansions on the target side.
    size_t view_size = nexp * sizeof(dcomplex_t);
    dcomplex_t *T1 = new dcomplex_t[nexp]();
    dcomplex_t *T2 = new dcomplex_t[nexp]();
    dcomplex_t *T3 = new dcomplex_t[nexp]();
    char *C1 = reinterpret_cast<char *>(T1);
    char *C2 = reinterpret_cast<char *>(T2);
    char *C3 = reinterpret_cast<char *>(T3);
    dcomplex_t *T[3] = {T1, T2, T3};
    char *C[3] = {C1, C2, C3};
    bool used[3] = {false, false, false};

    for (int i = 0; i < 3; ++i) {
      int tag = merge_and_shift_table[dx + 2][dy + 2][dz + 2][i];

      if (tag == -1) {
        break;
      }

      if (tag <= 1) {
        e2e(T[i], S_mz, dx, dy, 0, scale);
      } else if (tag <= 5) {
        e2e(T[i], S_my, dz, dx, 0, scale);
      } else if (tag <= 13) {
        e2e(T[i], S_mx, -dz, dy, 0, scale);
      } else if (tag <= 15) {
        e2e(T[i], S_pz, -dx, -dy, 0, scale);
      } else if (tag <= 19) {
        e2e(T[i], S_py, -dz, -dx, 0, scale);
      } else {
        e2e(T[i], S_px, dz, -dy, 0, scale);
      }

      views.add_view(tag, view_size, C[i]);
      used[i] = true;
    }

    if (used[1] == false) {
      delete [] T2;
    }

    if (used[2] == false) {
      delete [] T3;
    }

    expansion_t *retval = new expansion_t{views};
    return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> I_to_L(Index t_index, double t_size) const {
    // t_index and t_size is the index and size of the child
    // Compute child's center
    double h = t_size / 2;
    Point center = views_.center();
    double cx = center.x() + (t_index.x() % 2 == 0 ? -h : h);
    double cy = center.y() + (t_index.y() % 2 == 0 ? -h : h);
    double cz = center.z() + (t_index.z() % 2 == 0 ? -h : h);
    int to_child = 4 * (t_index.z() % 2) + 2 * (t_index.y() % 2) +
      (t_index.x() % 2);

    double scale = views_.scale() / 2;
    expansion_t *retval{new expansion_t{Point{cx, cy, cz},
                                        scale, kTargetPrimary}};

    int nexp = builtin_yukawa_table_->nexp(scale);

    dcomplex_t *E[28]{nullptr};
    for (int i = 0; i < 28; ++i) {
      E[i] = reinterpret_cast<dcomplex_t *>(views_.view_data(i));
    }
    dcomplex_t *L =
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    dcomplex_t *S = new dcomplex_t[nexp * 6]();
    dcomplex_t *S_mz = S;
    dcomplex_t *S_pz = S + nexp;
    dcomplex_t *S_my = S + 2 * nexp;
    dcomplex_t *S_py = S + 3 * nexp;
    dcomplex_t *S_mx = S + 4 * nexp;
    dcomplex_t *S_px = S + 5 * nexp;

    switch (to_child) {
    case 0:
      e2e(S_mz, E[uall], 0, 0, 3, scale);
      e2e(S_mz, E[u1234], 0, 0, 2, scale);
      e2e(S_pz, E[dall], 0, 0, 2, scale);

      e2e(S_my, E[nall], 0, 0, 3, scale);
      e2e(S_my, E[n1256], 0, 0, 2, scale);
      e2e(S_my, E[n12], 0, 0, 2, scale);
      e2e(S_py, E[sall], 0, 0, 2, scale);

      e2e(S_mx, E[eall], 0, 0, 3, scale);
      e2e(S_mx, E[e1357], 0, 0, 2, scale);
      e2e(S_mx, E[e13], 0, 0, 2, scale);
      e2e(S_mx, E[e1], 0, 0, 2, scale);
      e2e(S_px, E[wall], 0, 0, 2, scale);
      break;
    case 1:
      e2e(S_mz, E[uall], -1, 0, 3, scale);
      e2e(S_mz, E[u1234], -1, 0, 2, scale);
      e2e(S_pz, E[dall], 1, 0, 2, scale);

      e2e(S_my, E[nall], 0, -1, 3, scale);
      e2e(S_my, E[n1256], 0, -1, 2, scale);
      e2e(S_my, E[n12], 0, -1, 2, scale);
      e2e(S_py, E[sall], 0, 1, 2, scale);

      e2e(S_mx, E[eall], 0, 0, 2, scale);
      e2e(S_px, E[wall], 0, 0, 3, scale);
      e2e(S_px, E[w2468], 0, 0, 2, scale);
      e2e(S_px, E[w24], 0, 0, 2, scale);
      e2e(S_px, E[w2], 0, 0, 2, scale);
      break;
    case 2:
      e2e(S_mz, E[uall], 0, -1, 3, scale);
      e2e(S_mz, E[u1234], 0, -1, 2, scale);
      e2e(S_pz, E[dall], 0, 1, 2, scale);

      e2e(S_my, E[nall], 0, 0, 2, scale);
      e2e(S_py, E[sall], 0, 0, 3, scale);
      e2e(S_py, E[s3478], 0, 0, 2, scale);
      e2e(S_py, E[s34], 0, 0, 2, scale);

      e2e(S_mx, E[eall], 0, -1, 3, scale);
      e2e(S_mx, E[e1357], 0, -1, 2, scale);
      e2e(S_mx, E[e13], 0, -1, 2, scale);
      e2e(S_mx, E[e3], 0, -1, 2, scale);
      e2e(S_px, E[wall], 0, 1, 2, scale);
      break;
    case 3:
      e2e(S_mz, E[uall], -1, -1, 3, scale);
      e2e(S_mz, E[u1234], -1, -1, 2, scale);
      e2e(S_pz, E[dall], 1, 1, 2, scale);

      e2e(S_my, E[nall], 0, -1, 2, scale);
      e2e(S_py, E[sall], 0, 1, 3, scale);
      e2e(S_py, E[s3478], 0, 1, 2, scale);
      e2e(S_py, E[s34], 0, 1, 2, scale);

      e2e(S_mx, E[eall], 0, -1, 2, scale);
      e2e(S_px, E[wall], 0, 1, 3, scale);
      e2e(S_px, E[w2468], 0, 1, 2, scale);
      e2e(S_px, E[w24], 0, 1, 2, scale);
      e2e(S_px, E[w4], 0, 1, 2, scale);
      break;
    case 4:
      e2e(S_mz, E[uall], 0, 0, 2, scale);
      e2e(S_pz, E[dall], 0, 0, 3, scale);
      e2e(S_pz, E[d5678], 0, 0, 2, scale);

      e2e(S_my, E[nall], -1, 0, 3, scale);
      e2e(S_my, E[n1256], -1, 0, 2, scale);
      e2e(S_my, E[n56], -1, 0, 2, scale);
      e2e(S_py, E[sall], 1, 0, 2, scale);

      e2e(S_mx, E[eall], 1, 0, 3, scale);
      e2e(S_mx, E[e1357], 1, 0, 2, scale);
      e2e(S_mx, E[e57], 1, 0, 2, scale);
      e2e(S_mx, E[e5], 1, 0, 2, scale);
      e2e(S_px, E[wall], -1, 0, 2, scale);
      break;
    case 5:
      e2e(S_mz, E[uall], -1, 0, 2, scale);
      e2e(S_pz, E[dall], 1, 0, 3, scale);
      e2e(S_pz, E[d5678], 1, 0, 2, scale);

      e2e(S_my, E[nall], -1, -1, 3, scale);
      e2e(S_my, E[n1256], -1, -1, 2, scale);
      e2e(S_my, E[n56], -1, -1, 2, scale);
      e2e(S_py, E[sall], 1, 1, 2, scale);

      e2e(S_mx, E[eall], 1, 0, 2, scale);
      e2e(S_px, E[wall], -1, 0, 3, scale);
      e2e(S_px, E[w2468], -1, 0, 2, scale);
      e2e(S_px, E[w68], -1, 0, 2, scale);
      e2e(S_px, E[w6], -1, 0, 2, scale);
      break;
    case 6:
      e2e(S_mz, E[uall], 0, -1, 2, scale);
      e2e(S_pz, E[dall], 0, 1, 3, scale);
      e2e(S_pz, E[d5678], 0, 1, 2, scale);

      e2e(S_my, E[nall], -1, 0, 2, scale);
      e2e(S_py, E[sall], 1, 0, 3, scale);
      e2e(S_py, E[s3478], 1, 0, 2, scale);
      e2e(S_py, E[s78], 1, 0, 2, scale);

      e2e(S_mx, E[eall], 1, -1, 3, scale);
      e2e(S_mx, E[e1357], 1, -1, 2, scale);
      e2e(S_mx, E[e57], 1, -1, 2, scale);
      e2e(S_mx, E[e7], 1, -1, 2, scale);
      e2e(S_px, E[wall], -1, 1, 2, scale);
      break;
    case 7:
      e2e(S_mz, E[uall], -1, -1, 2, scale);
      e2e(S_pz, E[dall], 1, 1, 3, scale);
      e2e(S_pz, E[d5678], 1, 1, 2, scale);

      e2e(S_my, E[nall], -1, -1, 2, scale);
      e2e(S_py, E[sall], 1, 1, 3, scale);
      e2e(S_py, E[s3478], 1, 1, 2, scale);
      e2e(S_py, E[s78], 1, 1, 2, scale);

      e2e(S_mx, E[eall], 1, -1, 2, scale);
      e2e(S_px, E[wall], -1, 1, 3, scale);
      e2e(S_px, E[w2468], -1, 1, 2, scale);
      e2e(S_px, E[w68], -1, 1, 2, scale);
      e2e(S_px, E[w8], -1, 1, 2, scale);
      break;
    }

    e2l(S_mz, 'z', false, L);
    e2l(S_pz, 'z', true, L);
    e2l(S_my, 'y', false, L);
    e2l(S_py, 'y', true, L);
    e2l(S_mx, 'x', false, L);
    e2l(S_px, 'x', true, L);

    delete [] S;
    return std::unique_ptr<expansion_t>(retval);
  }

  void add_expansion(const expansion_t *temp1) {
    // This operation assumes that the views included in \p temp1 is a subset of
    // \p views_. No range checking performed.
    int count = temp1->views_.count();
    for (int i = 0; i < count; ++i) {
      int idx = temp1->views_.view_index(i);
      int size = temp1->views_.view_bytes(i) / sizeof(dcomplex_t);
      dcomplex_t *lhs = reinterpret_cast<dcomplex_t *>(views_.view_data(idx));
      dcomplex_t *rhs =
        reinterpret_cast<dcomplex_t *>(temp1->views_.view_data(i));

      for (int j = 0; j < size; ++j) {
        lhs[j] += rhs[j];
      }
    }
  }

  static void update_table(int n_digits, double domain_size,
                           const std::vector<double> &kernel_params) {
    update_yukawa_table(n_digits, domain_size, kernel_params[0]);
  }

  static void delete_table() { }

  static double compute_scale(Index index) {
    return builtin_yukawa_table_->scale(index.level());
  }

  static int weight_estimate(Operation op,
                             Index s = Index{}, Index t = Index{}) {
    int weight = 0;
    if (op == Operation::MtoI) {
      weight = 6;
    } else if (op == Operation::ItoI) {
      int weight = 0;
      int dx = s.x() - 2 * t.x();
      int dy = s.y() - 2 * t.y();
      int dz = s.z() - 2 * t.z();
      for (int i = 0; i < 3; ++i) {
        int tag = merge_and_shift_table[dx + 2][dy + 2][dz + 2][i];
        if (tag == -1) {
          break;
        }
        weight++;
      }
    } else {
      weight = 1;
    }
    return weight;
  }

private:
  ViewSet views_;

  void rotate_sph_z(const dcomplex_t *M, double alpha, dcomplex_t *MR) const {
    int p = builtin_yukawa_table_->p();
    // Compute exp(i * alpha)
    dcomplex_t ealpha = dcomplex_t{cos(alpha),  sin(alpha)};

    // Compute powers of exp(i * alpha)
    dcomplex_t *powers_ealpha = new dcomplex_t[p + 1];
    powers_ealpha[0] = dcomplex_t{1.0, 0.0};
    for (int j = 1; j <= p; ++j) {
      powers_ealpha[j] = powers_ealpha[j - 1] * ealpha;
    }

    int offset = 0;
    for (int n = 0; n <= p; ++n) {
      for (int m = 0; m <= n; ++m) {
        MR[offset] = M[offset] * powers_ealpha[m];
        offset++;
      }
    }

    delete [] powers_ealpha;
  }

  void rotate_sph_y(const dcomplex_t *M, const double *d,
                    dcomplex_t *MR) const {
    int p = builtin_yukawa_table_->p();
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

  void e2e(dcomplex_t *M, const dcomplex_t *W, int x, int y, int z,
           double scale) const {
    const dcomplex_t *xs = builtin_yukawa_table_->xs(scale);
    const dcomplex_t *ys = builtin_yukawa_table_->ys(scale);
    const double *zs = builtin_yukawa_table_->zs(scale);
    const int *m = builtin_yukawa_table_->m(scale);
    const int *sm = builtin_yukawa_table_->sm(scale);
    int s = builtin_yukawa_table_->s();

    int offset = 0;
    for (int k = 0; k < s; ++k) {
      double factor_z = zs[4 * k + z];
      for (int j = 0; j < m[k] / 2; ++j) {
        int idx = (sm[k] + j) * 7 + 3;
        dcomplex_t factor_x = xs[idx + x];
        dcomplex_t factor_y = ys[idx + y];
        M[offset] += W[offset] * factor_z * factor_y * factor_x;
        offset++;
      }
    }
  }

  void e2l(const dcomplex_t *E, char dir, bool sgn, dcomplex_t *L) const {
    // Note: this function is called on the parent node.
    double scale = views_.scale() / 2;;
    int p = builtin_yukawa_table_->p();
    int s = builtin_yukawa_table_->s();
    const int *m = builtin_yukawa_table_->m(scale);
    const int *sm = builtin_yukawa_table_->sm(scale);
    const int *smf = builtin_yukawa_table_->smf(scale);
    const int *f = builtin_yukawa_table_->f();
    const dcomplex_t *ealphaj = builtin_yukawa_table_->ealphaj(scale);
    double *legendre = new double[(p + 1) * (p + 2) / 2];
    const double *x = builtin_yukawa_table_->x();
    const double *w = builtin_yukawa_table_->w();
    double ld = builtin_yukawa_table_->lambda() *
      builtin_yukawa_table_->size(scale);
    const double *sqf = builtin_yukawa_table_->sqf();

    dcomplex_t *contrib = nullptr;
    dcomplex_t *W1 = new dcomplex_t[(p + 1) * (p + 2) / 2];
    dcomplex_t *W2 = new dcomplex_t[(p + 1) * (p + 2) / 2];

    for (int k = 0; k < s; ++k) {
      int mk2 = m[k] / 2;

      // Compute sum_{j=1}^m(k) W(k, j) e^{-i * m * alpha_j}
      dcomplex_t *z = new dcomplex_t[f[k] + 1];

      // m = 0
      z[0] = 0.0;
      for (int j = 1; j <= mk2; ++j) {
        int idx = sm[k] + j - 1;
        z[0] += (E[idx] + conj(E[idx]));
      }

      // m = 1, ..., f[k], where m is odd
      for (int m = 1; m <= f[k]; m += 2) {
        z[m] = 0.0;
        for (int j = 1; j <= mk2; ++j) {
          int widx = sm[k] + j - 1;
          int aidx = smf[k] + (j - 1) * f[k] + m - 1;
          z[m] += (E[widx] - conj(E[widx])) * conj(ealphaj[aidx]);
        }
      }

      // m = 2, ..., f[k], where m is even
      for (int m = 2; m <= f[k]; m += 2) {
        z[m] = 0.0;
        for (int j = 1; j <= mk2; ++j) {
          int widx = sm[k] + j - 1;
          int aidx = smf[k] + (j - 1) * f[k] + m - 1;
          z[m] += (E[widx] + conj(E[widx])) * conj(ealphaj[aidx]);
        }
      }

      legendre_Plm_gt1_scaled(p, 1 + x[k] / ld, scale, legendre);

      double factor = w[k] / m[k];
      for (int n = 0; n <= p; ++n) {
        int mmax = (n <= f[k] ? n : f[k]);
        for (int m = 0; m <= mmax; ++m) {
          W1[midx(n, m)] += z[m] * legendre[midx(n, m)] * factor;
        }
      }
    }

    // Scale the local expansion by
    // (2 * n + 1) * (n - m)! / (n + m)! * (-1)^n * i^m * pi / 2 / ld
    int offset = 0;
    double factor = M_PI_2 / ld;
    for (int n = 0; n <= p; ++n) {
      dcomplex_t power_I{1.0, 0.0};
      for (int m = 0; m <= n; ++m) {
        W1[offset++] *= sqf[midx(n, m)] * power_I * factor;
        power_I *= dcomplex_t{0.0, 1.0};
      }
      factor *= -1;
    }

    if (!sgn) {
      // If the exponential expansion is not along the positive axis
      // direction with respect to the source, flip the sign of the
      // converted L_n^m where n is odd
      int offset = 1;
      for (int n = 1; n <= p; n += 2) {
        for (int m = 0; m <= n; ++m) {
          W1[offset++] *= -1;
        }
        offset += (n + 2);
      }
    }

    if (dir == 'z') {
      contrib = W1;
    } else if (dir == 'y') {
      const double *d = builtin_yukawa_table_->dmat_plus(0.0);
      rotate_sph_y(W1, d, W2);
      rotate_sph_z(W2, M_PI / 2, W1);
      contrib = W1;
    } else if (dir == 'x') {
      const double *d = builtin_yukawa_table_->dmat_minus(0.0);
      rotate_sph_y(W1, d, W2);
      contrib = W2;
    }

    // Merge converted local expansion with the stored one
    offset = 0;
    for (int n = 0; n <= p; ++n) {
      for (int m = 0; m <= n; ++m) {
        L[offset] += contrib[offset];
        offset++;
      }
    }

    delete [] W1;
    delete [] W2;
    delete [] legendre;
  }
};


} // namespace dashmm

#endif // __DASHMM_YUKAWA_EXPANSION_H__




