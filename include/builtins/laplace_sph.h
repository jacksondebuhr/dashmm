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
#include <utility>
#include <initializer_list>

#include "dashmm/index.h"
#include "builtins/laplace_sph_table.h"
#include "dashmm/point.h"
#include "dashmm/types.h"
#include "dashmm/viewset.h"


namespace dashmm {

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

  LaplaceSPH(Point center, int n_digits, ExpansionRole role) 
    : views_{ViewSet{n_digits, role, center}} 
  {
    LaplaceSPHTableIterator entry = get_or_add_laplace_sph_table(n_digits);
    uLaplaceSPHTable &table = entry->second;

    // View size for each spherical harmonic expansion
    int p = table->p(); 
    int nsh = (p + 1) * (p + 2) / 2; 
    
    // View size for each exponential expansion 
    int nexp = table->nexp(); 

    if (role == kSourcePrimary || role == kTargetPrimary) {
      size_t bytes = sizeof(dcomplex_t) * nsh; 
      char *data = new char[bytes](); 
      views_.add_view(0, bytes, data); 
    } else if (role == kSourceIntermediate) {
      size_t bytes = sizeof(dcomplex_t) * nexp; 
      char *data = new char[bytes * 6](); 
      for (int i = 0; i < 6; ++i) 
        views_.add_view(i, bytes, data + bytes * i); 
    } else if (role == kTargetIntermediate) {
      size_t bytes = sizeof(dcomplex_t) * nexp; 
      char *data = new char[bytes * 28](); 
      for (int i = 0; i < 28; ++i) 
        views_.add_view(i, bytes, data + bytes * i);
    }
  }

  LaplaceSPH(const ViewSet &views) : views_{views} {
    int n_digits = views.n_digits(); 
    if (n_digits != -1) {
      get_or_add_laplace_sph_table(n_digits); 
    }
  }

  ~LaplaceSPH() {
    if (!views_.empty()) {
      for (auto it = views_.begin(); it != views_.end(); ++it) {
        delete [] (*it).data;
      }
    }
  }

  void release() {views_.clear();}

  bool valid(const ViewSet &view) const {
    // \p view is assumed to be a subset of \p views_ (no range checking
    // performed). The function returns true if and only if each entry in the
    // required subset is associated with some data. 
    bool is_valid = true; 
    for (auto it = view.begin(); it != view.end(); ++it) {
      int idx = (*it).index; 
      if (views_.view_data(idx) == nullptr) {
        is_valid = false; 
        break;
      }
    }
    return is_valid; 
  }

  int view_count() const { return views_.count(); }

  void get_views(ViewSet &view) const {
    // FIXME
    assert(view.count() < 2);
    if (view.count() > 0) {
      view.set_bytes(0, bytes_);
      view.set_data(0, (char *)data_);
    }
    view.set_n_digits(n_digits_);
    view.set_role(role_);
  }

  ViewSet get_all_views() const {
    // FIXME
    ViewSet retval{};
    retval.add_view(0);
    get_views(retval);
    return retval;
  }

  int accuracy() const {return views_.n_digits();}

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
                                      Source *last, double scale) const {
    int n_digits = views_.n_digits(); 
    expansion_t *retval{new expansion_t{center, n_digits, kSourcePrimary}};
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits);
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
          M[midx(n, m)] += q * powers_r[n] * powers_ephi[m] *
            legendre[midx(n, m)] * sqf[n - m] / sqf[n + m];
        }
      }
    }

    delete [] legendre;
    delete [] powers_r;
    delete [] powers_ephi;

    return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> S_to_L(Point center, Source *first,
                                      Source *last, double scale) const {
    int n_digits = views_.n_digits(); 
    expansion_t *retval{new expansion_t{center, n_digits, kTargetPrimary}};
    dcomplex_t *L = reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0)); 
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits);
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
          L[midx(n, m)] += q * powers_r[n] * powers_ephi[m] *
            legendre[midx(n, m)] * sqf[n - m] / sqf[n + m];
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
    Point center = views_.center(); 
    double px = center.x() + (from_child % 2 == 0 ? h : -h);
    double py = center.y() + (from_child % 4 <= 1 ? h : -h);
    double pz = center.z() + (from_child < 4 ? h : -h);

    int n_digits = views_.n_digits(); 
    expansion_t *retval{new expansion_t{Point{px, py, pz}, n_digits,
                                        kSourcePrimary}};

    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits);
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

    dcomplex_t *W1 = 
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    dcomplex_t *W2 = new dcomplex_t[(p + 1) * (p + 2) / 2];

    // Table of rotation angle about the z-axis, as an integer multiple of pi/4
    const int tab_alpha[8] = {1, 3, 7, 5, 1, 3, 7, 5};
    // Get rotation angle about the z-axis
    double alpha = tab_alpha[from_child] * M_PI_4;

    // Rotate the multipole expansion of the child box about z-axis
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
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

    int n_digits = views_.n_digits(); 
    expansion_t *retval{new expansion_t{Point{tx, ty, tz}, n_digits,
                                        kTargetPrimary}};
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits);
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
    dcomplex_t *W1 = 
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    dcomplex_t *W2 = new dcomplex_t[(p + 1) * (p + 2) / 2];

    // Compute the projection of t2s on the x-y plane
    const double proj = sqrt(t2s_x * t2s_x + t2s_y * t2s_y);

    // Handle of the multipole expansion
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(views_.view_data(0)); 

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
    Point center = views_.center(); 
    double h = t_size / 2;
    double cx = center.x() + (to_child % 2 == 0 ? -h : h);
    double cy = center.y() + (to_child % 4 <= 1 ? -h : h);
    double cz = center.z() + (to_child < 4 ? -h : h);

    int n_digits = views_.n_digits(); 
    expansion_t *retval{new expansion_t{Point{cx, cy, cz}, 
          n_digits, kTargetPrimary}};

    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits);
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
    
    dcomplex_t *W1 = 
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0)); 
    dcomplex_t *W2 = new dcomplex_t[(p + 1) * (p + 2) / 2];

    // Table of rotation angle about the z-axis as an integer multiple of pi / 4
    const int tab_alpha[8] = {1, 3, 7, 5, 1, 3, 7, 5};
    // Get rotation angle about the z-axis
    double alpha = tab_alpha[to_child] * M_PI_4;

    // Rotate the local expansion of the parent box about z-axis
    dcomplex_t *L = reinterpret_cast<dcomplex_t *>(views_.view_data(0)); 
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
    int n_digits = views_.n_digits(); 
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits);
    int p = table->p();
    const double *sqf = table->sqf();

    double *legendre = new double[(p + 1) * (p + 2) / 2];
    double *powers_r = new double[p + 1];
    dcomplex_t *powers_ephi = new dcomplex_t[p + 1];
    powers_ephi[0] = dcomplex_t{1.0, 0.0};
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(views_.view_data(0)); 

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
    int n_digits = views_.n_digits(); 
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits);
    int p = table->p();
    const double *sqf = table->sqf();

    double *legendre = new double[(p + 1) * (p + 2) / 2];
    double *powers_r = new double[p + 1];
    dcomplex_t *powers_ephi = new dcomplex_t[p + 1];
    powers_r[0] = 1.0;
    powers_ephi[0] = dcomplex_t{1.0, 0.0};

    dcomplex_t *L = reinterpret_cast<dcomplex_t *>(views_.view_data(0)); 

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

  std::unique_ptr<expansion_t> M_to_I(Index s_index) const {   
    int n_digits = views_.n_digits(); 
    expansion_t *retval{new expansion_t{views_.center(), 
          n_digits, kSourceIntermediate}}; 
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(views_.view_data(0)); 
    dcomplex_t *E = 
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0)); 

    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits); 
    int p = table->p(); 
    int s = table->s(); 
    int nsh = (p + 1) * (p + 2) / 2; 
    int nexp = table->nexp(); 
    const double *weight = table->weight(); 
    const int *m_ = table->m(); 
    const int *f_ = table->f(); 
    const int *smf_ = table->smf(); 
    const double *d1 = table->dmat_plus(0.0); 
    const double *d2 = table->dmat_minus(0.0); 
    const double *lambdaknm = table->lambdaknm(); 
    const dcomplex_t *ealphaj = table->ealphaj(); 

    // Allocate scratch space to handle x- and y-direction exponential expansions
    dcomplex_t *W1 = new dcomplex_t[(p + 1) * (p + 2) / 2]; 
    dcomplex_t *W2 = new dcomplex_t[(p + 1) * (p + 2) / 2]; 

    // Addresses of exponential expansions in the positive axis direction
    dcomplex_t *EP[3] = {E, E + 2 * nexp, E + 4 * nexp}; 
    // Addresses of exponential expansions in the negative axis direction
    dcomplex_t *EM[3] = {E + nexp, E + 3 * nexp, E + 5 * nexp}; 

    // Setup y-direction. Rotate the multipole expansion M about z-axis by -pi /
    // 2, making (x, y, z) frame (-y, x, z). Next, rotate it again about the new
    // y axis by -pi / 2. The (-y, x, z) in the first rotated frame becomes (z,
    // x, y) in the final frame. 
    retval->rotate_sph_z(M, -M_PI / 2, W1); 
    retval->rotate_sph_y(W1, d2, W2); 
    
    // Setup x-direction. Rotate the multipole expansion M about y axis by pi /
    // 2. This makes (x, y, z) frame into (-z, y, x). 
    retval->rotate_sph_y(M, d1, W1); 

    // Addresses of the spherical harmonic expansions 
    const dcomplex_t *SH[3] = {&W1[0], &W2[0], M}; 

    for (int dir = 0; dir <= 2; ++dir) {
      int offset = 0; 
      for (int k = 0; k < s ; ++k) {
        double weight = weight[k] / m_[k]; 
                
        // Compute sum_{n = m}^p M_n^m * lambda_k^n / sqrt((n+m)! * (n - m)!)
        // z1 handles M_n^m where n is even 
        dcomplex_t *z1 = new dcomplex_t[f_[k] + 1]; 
        // z2 handles M_n^m where n is odd 
        dcomplex_t *z2 = new dcomplex_t[f_[k] + 1]; 
        
        // Process M_n^0 terms
        z1[0] = 0; 
        z2[0] = 0; 
        for (int n = 0; n <= p; n += 2) 
          z1[0] += SH[dir][midx(n, 0)] * lambdaknm[k * nsh + midx(n, 0)]; 
        for (int n = 1; n <= p; n += 2) 
          z2[0] += SH[dir][midx(n, 0)] * lambdaknm[k * nsh + midx(n, 0)]; 
      
        // Process M_n^m terms for nonzero m
        for (int m = 1; m <= f_[k]; m += 2) {
          z1[m] = 0; 
          z2[m] = 0; 
          for (int n = m; n <= p; n += 2) 
            z2[m] += SH[dir][midx(n, m)] * lambdaknm[k * nsh + midx(n, m)]; 
          for (int n = m + 1; n <= p_; n += 2) 
            z1[m] += SH[dir][midx(n, m)] * lambdaknm[k * nsh + midx(n, m)]; 
        }

        for (int m = 2; m <= f_[k]; m += 2) {
          z1[m] = 0;
          z2[m] = 0;
          for (int n = m; n <= p; n += 2) 
            z1[m] += SH[dir][midx(n, m)] * lambdaknm[k * nsh + midx(n, m)]; 
          for (int n = m + 1; n <= p_; n += 2) 
            z2[m] += SH[dir][midx(n, m)] * lambdaknm[k * nsh + midx(n, m)]; 
        }
        
        // Compute W(k, j) 
        for (int j = 1; j <= m_[k] / 2; ++j) {
          dcomplex_t up = z1[0] + z2[0]; // accumulate +dir
          dcomplex_t dn = z1[0] - z2[0]; // accumulate -dir
          dcomplex_t power_I {0.0, 1.0}; 
          for (int m = 1; m <= f_[k]; ++m) {
            int idx = smf_[k] + (j - 1) * f_[k] + (m - 1); 
            up += 2 * real(ealphaj[idx]) * (z1[m] + z2[m]) * power_I; 
            dn += 2 * real(ealphaj[idx]) * (z1[m] - z2[m]) * power_I; 
            power_I *= dcomplex_t{0.0, 1.0};
          }
          EP[dir][offset] = weight * up;
          EM[dir][offset] = weight * dn;
          offset++; 
        }
        
        delete [] z1; 
        delete [] z2; 
      }
    }

    delete [] W1; 
    delete [] W2; 

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
    int n_digits = views_.n_digits(); 
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits); 
    int nexp = table->nexp(); 
    const dcomplex_t *S = 
      reinterpret_cast<dcomplex_t *>(views_.view_data(0)); 
    const dcomplex_t *S_px = S; // +x direction
    const dcomplex_t *S_mx = S + nexp; // -x direction
    const dcomplex_t *S_py = S + nexp * 2; // +y direction
    const dcomplex_t *S_my = S + nexp * 3; // -y direction
    const dcomplex_t *S_pz = S + nexp * 4; // +z direction
    const dcomplex_t *S_mz = S + nexp * 5; // -z direction

    ViewSet views{n_digits, kTargetIntermediate, Point{px, py, pz}}; 

    // Each S is at most going to contribute 3 views to the exponential
    // expansions on the target side. 
    dcomplex_t *work = new dcomplex_t[nexp * 3](); 
    dcomplex_t *T1 = work; 
    dcomplex_t *T2 = work + nexp; 
    dcomplex_t *T3 = work + nexp * 2; 

    // Produce views needed on the target side 
    if (dz == 3) {
      this->e2e(T1, {std::make_pair(S_mz, Index{dx, dy, 0})}); 
      views.add_view(uall, nexp, T1); 
    } else if (dz == -2) {
      this->e2e(T1, {std::make_pair(S_pz, Index{-dx, -dy, 0})}); 
      views.add_view(dall, nexp, T1); 
    } else if (dy == 3) {
      this->e2e(T1, {std::make_pair(S_my, Index{dz, dx, 0})}); 
      views.add_view(nall, nexp, T1);
    } else if (dy == -2) {
      this->e2e(T1, {std::make_pair(S_py, Index{-dz, -dx, 0})}); 
      views.add_view(sall, nexp, T1); 
    } else if (dx == 3) {
      this->e2e(T1, {std::make_pair(S_mx, Index{-dz, dy, 0})}); 
      views.add_view(eall, nexp, T1); 
    } else if (dx == -2) {
      this->e2e(T1, {std::make_pair(S_px, Index{dz, -dy, 0})}); 
      views.add_view(wall, nexp, T1); 
    } else {
      if (dz == 2) {
        this->e2e(T1, {std::make_pair(S_mz, Index{dx, dy, 0})});
        views.add_view(u1234, nexp, T1); 

        if (dy == -1) {
          this->e2e(T2, {std::make_pair(S_py, Index{-dz, -dx, 0})}); 
          views.add_view(s78, nexp, T2); 
        
          if (dx == -1) {
            this->e2e(T3, {std::make_pair(S_px, Index{dz, -dy, 0})}); 
            views.add_view(w6, nexp, T3); 
          } else if (dx == 2) {
            this->e2e(T3, {std::make_pair(S_mx, Index{-dz, dy, 0})}); 
            views.add_view(e5, nexp, T3); 
          }
        } else if (dy == 2) {
          this->e2e(T2, {std::make_pair(S_my, Index{dz, dx, 0})});
          views.add_view(n56, nexp, T2); 
        
          if (dx == -1) {
            this->e2e(T3, {std::make_pair(S_px, Index{dz, -dy, 0})});
            views.add_view(w8, nexp, T3); 
          } else if (dx == 2) {
            this->e2e(T3, {std::make_pair(S_mx, Index{-dz, dy, 0})});
            views.add_view(e7, nexp, T3); 
          } 
        } else {
          if (dx == -1) {            
            this->e2e(T2, {std::make_pair(S_px, Index{dz, -dy, 0})});
            views.add_view(w68, nexp, T2); 
          } else if (dx == 2) {
            this->e2e(T2, {std::make_pair(S_mx, Index{-dz, dy, 0})});
            views.add_view(e57, nexp, T2); 
          }
        }
      } else if (dz == -1) {
        this->e2e(T1, {std::make_pair(S_pz, Index{-dx, -dy, 0})});
        views.add_view(d5678, nexp, T1);

        if (dy == -1) {
          this->e2e(T2, {std::make_pair(S_py, Index{-dz, -dx, 0})});
          views.add_view(s34, nexp, T2); 

          if (dx == -1) {
            this->e2e(T3, {std::make_pair(S_px, Index{dz, -dy, 0})});
            views.add_view(w2, nexp, T3); 
          } else if (dx == 2) {
            this->e2e(T3, {std::make_pair(S_mx, Index{-dz, dy, 0})});
            views.add_view(e1, nexp, T3); 
          } 
        } else if (dy == 2) {
          this->e2e(T2, {std::make_pair(S_my, Index{dz, dx, 0})});
          views.add_view(n12, nexp, T2); 

          if (dx == -1) {
            this->e2e(T3, {std::make_pair(S_px, Index{dz, -dy, 0})});
            views.add_view(w4, nexp, T3); 
          } else if (dx == 2) {
            this->e2e(T3, {std::make_pair(S_mx, Index{-dz, dy, 0})});
            views.add_view(e3, nexp, T3); 
          } 
        } else {
          if (dx == -1) {
            this->e2e(T2, {std::make_pair(S_px, Index{dz, -dy, 0})});
            views.add_view(w24, nexp, T2); 
          } else if (dx == 2) {
            this->e2e(T2, {std::make_pair(S_mx, Index{-dz, dy, 0})});
            views.add_view(e13, nexp, T2); 
          }
        }
      } else { 
        if (dy == -1) {
          this->e2e(T1, {std::make_pair(S_py, Index{-dz, -dx, 0})});
          views.add_view(s3478, nexp, T1); 

          if (dx == -1) {
            this->e2e(T2, {std::make_pair(S_px, Index{dz, -dy, 0})});
            this->e2e(T3, {std::make_pair(S_px, Index{dz, -dy, 0})});
            views.add_view(w2, nexp, T2); 
            views.add_view(w6, nexp, T3); 
          } else if (dx == 2) {
            this->e2e(T2, {std::make_pair(S_mx, Index{-dz, dy, 0})});
            this->e2e(T3, {std::make_pair(S_mx, Index{-dz, dy, 0})});
            views.add_view(e1, nexp, T2); 
            views.add_view(e5, nexp, T3); 
          } 
        } else if (dy == 2) {
          this->e2e(T1, {std::make_pair(S_my, Index{dz, dx, 0})});
          views.add_view(n1256, nexp, T1); 

          if (dx == -1) {
            this->e2e(T2, {std::make_pair(S_px, Index{dz, -dy, 0})});
            this->e2e(T3, {std::make_pair(S_px, Index{dz, -dy, 0})});
            views.add_view(w4, nexp, T2); 
            views.add_view(w8, nexp, T3); 
          } else if (dx == 2) {
            this->e2e(T2, {std::make_pair(S_mx, Index{-dz, dy, 0})});
            this->e2e(T3, {std::make_pair(S_mx, Index{-dz, dy, 0})});
            views.add_view(e3, nexp, T2); 
            views.add_view(e7, nexp, T3); 
          }
        } else { 
          if (dx == -1) {
            this->e2e(T2, {std::make_pair(S_px, Index{dz, -dy, 0})});
            views.add_view(w2468, nexp, T2); 
          } else if (dx == 2) {
            this->e2e(T2, {std::make_pair(S_mx, Index{-dz, dy, 0})});
            views.add_view(e1357, nexp, T2); 
          }
        }
      }
    }
  
    expansion_t retval{views};       
    return std::unique_ptr<expansion_t>{&retval}; 
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

    int n_digits = views_.n_digits(); 
    expansion_t *retval{new expansion_t{Point{cx, cy, cz}, 
          n_digits_, kTargetPrimary}}; 
    
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits);
    int nexp = table->nexp(); 

    dcomplex_t *E = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    dcomplex_t *L = 
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0)); 
    dcomplex_t *S = new dcomplex[nexp * 6](); 
    dcomplex_t *S_mz = S; 
    dcomplex_t *S_pz = S + nexp; 
    dcomplex_t *S_my = S + 2 * nexp; 
    dcomplex_t *S_py = S + 3 * nexp; 
    dcomplex_t *S_mx = S + 4 * nexp; 
    dcomplex_t *S_px = S + 5 * nexp; 

    int to_child = (t_index.x() % 2) + 2 * (t_index.y() % 2) + 
      4 * (t_index.z() % 2); 

    switch (to_child) {
    case 0: 
      this->e2e(S_mz, {std::make_pair(&M[uall * nexp], Index{0, 0, 3}), 
            std::make_pair(&M[u1234 * nexp], Index{0, 0, 2})});
      this->e2e(S_pz, {std::make_pair(&M[dall * nexp], Index{0, 0, 2})});
      
      this->e2e(S_my, {std::make_pair(&M[nall * nexp], Index{0, 0, 3}),
            std::make_pair(&M[n1256 * nexp], Index{0, 0, 2}), 
            std::make_pair(&M[n12 * nexp], Index{0, 0, 2})});
      this->e2e(S_py, {std::make_pair(&M[sall * nexp], Index{0, 0, 2})});
      
      this->e2e(S_mx, {std::make_pair(&M[eall * nexp], Index{0, 0, 3}), 
            std::make_pair(&M[e1357 * nexp], Index{0, 0, 2}), 
            std::make_pair(&M[e13 * nexp], Index{0, 0, 2}), 
            std::make_pair(&M[e1 * nexp], Index{0, 0, 2})});
      this->e2e(S_px, {std::make_pair(&M[wall * nexp], Index{0, 0, 2})});
      break; 
    case 1: 
      this->e2e(S_mz, {std::make_pair(&M[uall * nexp], Index{-1, 0, 3}), 
            std::make_pair(&M[u1234 * nexp], Index{-1, 0, 2})}); 
      this->e2e(S_pz, {std::make_pair(&M[dall * nexp], Index{1, 0, 2})}); 
      
      this->e2e(S_my, {std::make_pair(&M[nall * nexp], Index{0, -1, 3}), 
            std::make_pair(&M[n1256 * nexp], Index{0, -1, 2}), 
            std::make_pair(&M[n12 * nexp], Index{0, -1, 2})}); 
      this->e2e(S_py, {std::make_pair(&M[sall * nexp], Index{0, 1, 2})});
      
      this->e2e(S_mx, {std::make_pair(&M[eall * nexp], Index{0, 0, 2})});
      this->e2e(S_px, {std::make_pair(&M[wall * nexp], Index{0, 0, 3}), 
            std::make_pair(&M[w2468 * nexp], Index{0, 0, 2}), 
            std::make_pair(&M[w24 * nexp], Index{0, 0, 2}), 
            std::make_pair(&M[w2 * nexp], Index{0, 0, 2})});
      break;
    case 2: 
      this->e2e(S_mz, {std::make_pair(&M[uall * nexp], Index{0, -1, 3}), 
            std::make_pair(&M[u1234 * nexp], Index{0, -1, 2})});
      this->e2e(S_pz, {std::make_pair(&M[dall * nexp], Index{0, 1, 2})});
      
      this->e2e(S_my, {std::make_pair(&M[nall * nexp], Index{0, 0, 2})});
      this->e2e(S_py, {std::make_pair(&M[sall * nexp], Index{0, 0, 3}), 
            std::make_pair(&M[s3478 * nexp], Index{0, 0, 2}), 
            std::make_pair(&M[s34 * nexp], Index{0, 0, 2})}); 
      
      this->e2e(S_mx, {std::make_pair(&M[eall * nexp], Index{0, -1, 3}), 
            std::make_pair(&M[e1357 * nexp], Index{0, -1, 2}), 
            std::make_pair(&M[e13 * nexp], Index{0, -1, 2}), 
            std::make_pair(&M[e3 * nexp], Index{0, -1, 2})}); 
      this->e2e(S_px, {std::make_pair(&M[wall * nexp], Index{0, 1, 2})}); 
      break; 
    case 3: 
      this->e2e(S_mz, {std::make_pair(&M[uall * nexp], Index{-1, -1, 3}), 
            std::make_pair(&M[u1234 * nexp], Index{-1, -1, 2})});
      this->e2e(S_pz, {std::make_pair(&M[dall * nexp], Index{1, 1, 2})}); 
      
      this->e2e(S_my, {std::make_pair(&M[nall * nexp], Index{0, -1, 2})}); 
      this->e2e(S_py, {std::make_pair(&M[sall * nexp], Index{0, 1, 3}), 
            std::make_pair(&M[s3478 * nexp], Index{0, 1, 2}), 
            std::make_pair(&M[s34 * nexp], Index{0, 1, 2})}); 

      this->e2e(S_mx, {std::make_pair(&M[eall * nexp], Index{0, -1, 2})}); 
      this->e2e(S_px, {std::make_pair(&M[wall * nexp], Index{0, 1, 3}), 
            std::make_pair(&M[w2468 * nexp], Index{0, 1, 2}), 
            std::make_pair(&M[w24 * nexp], Index{0, 1, 2}), 
            std::make_pair(&M[w4 * nexp], Index{0, 1, 2})}); 
      break; 
    case 4: 
      this->e2e(S_mz, {std::make_pair(&M[uall * nexp], Index{0, 0, 2})}); 
      this->e2e(S_pz, {std::make_pair(&M[dall * nexp], Index{0, 0, 3}), 
            std::make_pair(&M[d5678 * nexp], Index{0, 0, 2})}); 
      
      this->e2e(S_my, {std::make_pair(&M[nall * nexp], Index{-1, 0, 3}), 
            std::make_pair(&M[n1256 * nexp], Index{-1, 0, 2}), 
            std::make_pair(&M[n56 * nexp], Index{-1, 0, 2})}); 
      this->e2e(S_py, {std::make_pair(&M[sall * nexp], Index{1, 0, 2})}); 
      
      this->e2e(S_mx, {std::make_pair(&M[eall * nexp], Index{1, 0, 3}), 
            std::make_pair(&M[e1357 * nexp], Index{1, 0, 2}), 
            std::make_pair(&M[e57 * nexp], Index{1, 0, 2}), 
            std::make_pair(&M[e5 * nexp], Index{1, 0, 2})}); 
      this->e2e(S_px, {std::make_pair(&M[wall * nexp], Index{-1, 0, 2})}); 
      break; 
    case 5: 
      this->e2e(S_mz, {std::make_pair(&M[uall * nexp], Index{-1, 0, 2})}); 
      this->e2e(S_pz, {std::make_pair(&M[dall * nexp], Index{1, 0, 3}), 
            std::make_pair(&M[d5678 * nexp], Index{1, 0, 2})}); 
      
      this->e2e(S_my, {std::make_pair(&M[nall * nexp], Index{-1, -1, 3}), 
            std::make_pair(&M[n1256 * nexp], Index{-1, -1, 2}), 
            std::make_pair(&M[n56 * nexp], Index{-1, -1, 2})}); 
      this->e2e(S_py, {std::make_pair(&M[sall * nexp], Index{1, 1, 2})}); 
      
      this->e2e(S_mx, {std::make_pair(&M[eall * nexp], Index{1, 0, 2})});  
      this->e2e(S_px, {std::make_pair(&M[wall * nexp], Index{-1, 0, 3}), 
            std::make_pair(&M[w2468 * nexp], Index{-1, 0, 2}), 
            std::make_pair(&M[w68 * nexp], Index{-1, 0, 2}), 
            std::make_pair(&M[w6 * nexp], Index{-1, 0, 2})}); 
      break;
    case 6: 
      this->e2e(S_mz, {std::make_pair(&M[uall * nexp], Index{0, -1, 2})});
      this->e2e(S_pz, {std::make_pair(&M[dall * nexp], Index{0, 1, 3}), 
            std::make_pair(&M[d5678 * nexp], Index{0, 1, 2})}); 
      
      this->e2e(S_my, {std::make_pair(&M[nall * nexp], Index{-1, 0, 2})});  
      this->e2e(S_py, {std::make_pair(&M[sall * nexp], Index{1, 0, 3}), 
            std::make_pair(&M[s3478 * nexp], Index{1, 0, 2}), 
            std::make_pair(&M[s78 * nexp], Index{1, 0, 2})}); 
      
      this->e2e(S_mx, {std::make_pair(&M[eall * nexp], Index{1, -1, 3}), 
            std::make_pair(&M[e1357 * nexp], Index{1, -1, 2}), 
            std::make_pair(&M[e57 * nexp], Index{1, -1, 2}), 
            std::make_pair(&M[e7 * nexp], Index{1, -1, 2})}); 
      this->e2e(S_px, {std::make_pair(&M[wall * nexp], Index{-1, 1, 2})}); 
      break; 
    case 7: 
      this->e2e(S_mz, {std::make_pair(&M[uall * nexp], Index{-1, -1, 2})}); 
      this->e2e(S_pz, {std::make_pair(&M[dall * nexp], Index{1, 1, 3}), 
            std::make_pair(&M[d5678 * nexp], Index{1, 1, 2})}); 
      
      this->e2e(S_my, {std::make_pair(&M[nall * nexp], Index{-1, -1, 2})});  
      this->e2e(S_py, {std::make_pair(&M[sall * nexp], Index{1, 1, 3}), 
            std::make_pair(&M[s3478 * nexp], Index{1, 1, 2}), 
            std::make_pair(&M[s78 * nexp], Index{1, 1, 2})}); 
      
      this->e2e(S_mx, {std::make_pair(&M[eall * nexp], Index{1, -1, 2})}); 
      this->e2e(S_px, {std::make_pair(&M[wall * nexp], Index{-1, 1, 3}), 
            std::make_pair(&M[w2468 * nexp], Index{-1, 1, 2}), 
            std::make_pair(&M[w68 * nexp], Index{-1, 1, 2}), 
            std::make_pair(&M[w8 * nexp], Index{-1, 1, 2})}); 
      break;
    }

    this->e2l(S_mz, 'z', false, L); 
    this->e2l(S_pz, 'z', true, L); 
    this->e2l(S_my, 'y', false, L); 
    this->e2l(S_py, 'y', true, L); 
    this->e2l(S_mx, 'x', false, L); 
    this->e2l(S_px, 'x', true, L); 

    return std::unique_ptr<expansion_t>(retval); 
  }

  void add_expansion(const expansion_t *temp1) {
    // This operation assumes that the views included in \p temp1 is a subset of
    // \p views_. No range checking performed. 
    int count = temp1->views_.count(); 
    for (int i = 0; i < count; ++i) {
      int idx = temp1->views_.view_index(i); 
      int size = temp1->views_.view_size(i); 
      dcomplex_t *lhs = reinterpret_cast<dcomplex_t *>(views_.view_data(idx));
      dcomplex_t *rhs = 
        reinterpret_cast<dcomplex_t *>(temp1->views_.view_data(i));

      for (int j = 0; j < size; ++j) 
        lhs[j] += rhs[j];
    }
  }

 private:
  ViewSet views_; 

  void rotate_sph_z(const dcomplex_t *M, double alpha, dcomplex_t *MR) {
    int n_digits = views_.n_digits(); 
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits);
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
    int n_digits = views_.n_digits(); 
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits);
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
    int n_digits = views_.n_digits(); 
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits);
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
    int n_digits = views_.n_digits(); 
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits);
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

  void e2e(dcomplex_t *M, std::initializer_list<e2e_t> inputs) {
    int n_digits = views_.n_digits(); 
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits); 
    const dcomplex_t *xs = table->xs(); 
    const dcomplex_t *ys = table->ys(); 
    const double *zs = table->zs(); 
    int s = table->s(); 
    const int *m = table->m(); 
    const int *sm = table->sm(); 

    for (auto it = inputs.begin(); it != inputs.end(); ++it) {
      const dcomplex_t *W = std::get<0>(*it); 
      Index I = std::get<1>(*it); 
      
      int offset = 0; 
      for (int k = 0; k < s; ++k) {
      // Shifting factor in z direction
        double factor_z = zs[4 * k + I.z()]; 
        for (int j = 0; j < m[k] / 2; ++j) {
          int sidx = (sm[j] + j) * 7 + 3; 
          // Shifting factor in x direction
          dcomplex_t factor_x = xs[sidx + I.x()]; 
          // Shifting factor in y direction
          dcomplex_t factor_y = ys[sidx + I.y()]; 
          M[offset] += W[offset] * factor_z * factor_y * factor_x; 
          offset++; 
        }
      }
    }
  }

  void e2l(const dcomplex_t *E, char dir, bool sgn, dcomplex_t *L) {
    int n_digits = views_.n_digits(); 
    uLaplaceSPHTable &table = builtin_laplace_table_.at(n_digits); 
    const double *sqf = table->sqf(); 
    const dcomplex_t *ealphaj = table->ealphaj(); 
    const int *m_ = table->m(); 
    const int *sm_ = table->sm(); 
    const int *f_ = table->f(); 
    const int *smf_ = table->smf(); 
    int p = table->p(); 
    int s = table->s(); 

    dcomplex_t *contrib = nullptr; 
    dcomplex_t *W1 = new dcomplex_t[(p + 1) * (p + 2) / 2]; 
    dcomplex_t *W2 = new dcomplex_t[(p + 1) * (p + 2) / 2];

    for (int k = 0; k < s; ++k) {
      int Mk2 = m_[k] / 2; 

      // Compute sum_{j = 1}^{M(k)} W(k, j) exp(-i * m * alpha_j) 
      dcomplex_t *z = new dcomplex_t[f_[k] + 1]; 

      // m = 0 
      z[0] = 0.0; 
      for (int j = 1; j <= Mk2; ++j) {
        int idx = sm_[k] + j - 1; 
        z[0] += (E[idx] + conj(E[idx]));
      }

      // m = 1, ..., F(k), where m is odd 
      for (int m = 1; m <= f_[k]; m += 2) {
        z[m] = 0.0; 
        for (int j = 1; j <= Mk2; ++j) {
          int idx1 = sm_[k] + j - 1; 
          int idx2 = smf_[k] + (j - 1) * f_[k] + m - 1; 
          z[m] += (E[idx1] - conj(E[idx1])) * conj(ealphaj[idx2]);          
        }
      }
        
      // m = 2, ..., F(k), where m is even 
      for (int m = 2; m <= f_[k]; m += 2) {
        z[m] = 0.0; 
        for (int j = 1; j <= Mk2; ++j) {
          int idx1 = sm_[k] + j - 1; 
          int idx2 = smf_[k] + (j - 1) * f_[k] + m - 1; 
          z[m] += (E[idx1] + conj(E[idx1])) * conj(ealphaj[idx2]);
        }
      }

      // Compute lambda_k's contribution 
      double power_lambdak = 1.0; // (-lambda_k)^n 
      for (int n = 0; n <= p; ++n) {
        int mmax = fmin(n, f_[k]); 
        for (int m = 0; m <= mmax; ++m) 
          W1[midx(n, m)] += power_lambdak * z[m]; 
        power_lambdak *= -lambda_[k];
      }
      delete [] z; 
    }

    if (!sgn) {
      // If the exponential expansion is not along the positive direction of the
      // axis with respect to the source, flip the sign of the converted L_n^m
      // terms where n is odd. 
      int offset = 1; // address of L_1^0 
      for (int n = 1; n <= p; n += 2) {
        for (int m = 0; m <= n; ++m) {
          W1[offset++] *= -1; 
        }
        offset += (n + 2); // skip the storage for the even value of n 
      }
    }

    // Scale the local expansion by i^m / sqrt((n - m)!(n + m)!)
    int offset = 0; 
    for (int n = 0; n <= p; ++n) {
      dcomplex_t power_I{1.0, 0.0}; 
      for (int m = 0; m <= n; ++m) {
        W1[offset++] *= power_I / sqf[n - m] / sqf[n + m]; 
        power_I *= dcomplex_t{0.0, 1.0}; 
      }
    }

    if (dir == 'z') {
      contrib = W1; 
    } else if (dir == 'y') {
      const double *d = table->dmat_plus(0.0); 
      this->rotate_sph_y(W1, d, W2); 
      this->rotate_sph_z(W2, M_PI / 2, W1); 
      contrib = W1; 
    } else if (dir == 'x') {
      const double *d = table->dmat_minus(0.0); 
      this->rotate_sph_y(W1, d, W2); 
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
  }
};


} // namespace dashmm

#endif // __DASHMM_LAPLACE_SPH_EXPANSION_H__
