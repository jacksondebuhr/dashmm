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

#ifndef __DASHMM_HELMHOLTZ_EXPANSION_H__
#define __DASHMM_HELMHOLTZ_EXPANSION_H__

/// \file 
/// \brief Declaration of Helmholtz (low-frequency) 

#include <cassert>
#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <vector>

#include "dashmm/index.h"
#include "builtins/helmholtz_table.h"
#include "dashmm/point.h"
#include "dashmm/types.h"
#include "dashmm/viewset.h"

namespace dashmm {

/// Helmholtz kernel Spherical Harmonic expansion 
///
/// This expansion is of the Helmholtz kernel about the center of the node
/// containing the represented sources. 
/// 
/// this class is a template with parameters for the source and target types. 
///
/// Source must define a double valued 'charge' member to be used with
/// Helmholtz. Target must define a std::complex<double> valued 'phi' member to
/// be used with Helmholtz. 
template <typename Source, typename Target>
class Helmholtz {
public: 
  using source_t = Source; 
  using target_t = Target; 
  using expansion_t = Helmholtz<Source, Target>; 

  Helmholtz(Point center, double scale, ExpansionRole role) 
    : views_{ViewSet{role, center, scale}} {
    
    // View size for each spherical harmonic expansion
    int p = builtin_helmholtz_table_->p(); 
    int n_e = builtin_helmholtz_table_->n_e(scale); 
    int n_p = builtin_helmholtz_table_->n_p(scale); 

    if (role == kSourcePrimary) {
      size_t bytes = sizeof(dcomplex_t) * (p + 1) * (p + 2) / 2; 
      char *data = new char[bytes](); 
      views_.add_view(0, bytes, data); 
    } else if (role == kTargetPrimary) {
      size_t bytes = sizeof(dcomplex_t) * (p + 1) * (p + 1); 
      char *data = new char[bytes](); 
      views_.add_view(0, bytes, data); 
    } else if (role == kSourceIntermediate) {
      // On the source side, propagating waves along positive and negative
      // directions of a given axis are conjugate of each other. As a result,
      // only the one along the positive direction is saved. 
      size_t bytes_p = sizeof(dcomplex_t) * n_p; 
      size_t bytes_e = sizeof(dcomplex_t) * n_e; 

      // On the source side, one have
      // (1) Prop(k, j) = conj(Prop(k, j + m_p(k) / 2) 
      // (2) Evan(k, j) = conj(Evan(k, j + m_e(k) / 2)
      // Only half of the coefficients are saved for this reason.
      for (int i = 0; i < 3; ++i) {
        int j = 3 * i; 

        // Propagating wave 
        char *data1 = new char[bytes_p](); 
        views_.add_view(j, bytes_p, data1); 

        // Evanescent wave positive axis
        char *data2 = new char[bytes_e](); 
        views_.add_view(j + 1, bytes_e, data2); 

        // Evanescent wave negative axis 
        char *data3 = new char[bytes_e](); 
        views_.add_view(j + 2, bytes_e, data3); 
      }
    } else if (role == kTargetIntermediate) {
      // On the target side, propgating wave loses symmetry, and one needs to
      // store Prop(k, j) for 1 <= j <= m_p(k). Evanescent wave still maintains
      // the symmetry, so only half of the Evan(k, j) are saved. 
      size_t bytes_p = sizeof(dcomplex_t) * n_p * 2; 
      size_t bytes_e = sizeof(dcomplex_t) * n_e; 
      
      for (int i = 0; i < 28; ++i) {
        int j = 2 * i; 
        
        // Propagating wave 
        char *data1 = new char[bytes_p](); 
        views_.add_view(j, bytes_p, data1); 

        // Evanescent wave 
        char *data2 = new char[bytes_e](); 
        views_.add_view(j + 1, bytes_e, data2); 
      }
    }
  }

  Helmholtz(const ViewSet &views) : views_{views} { }

  ~Helmholtz() {
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
    int p = builtin_helmholtz_table_->p(); 
    const double *sqf = builtin_helmholtz_table_->sqf(); 
    double omega = builtin_helmholtz_table_->omega(); 

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
                         dcomplex_t{dist.x() / proj, -dist.y() / proj}); 

      // Compute powers of exp(-i * phi)
      powers_ephi[0] = 1.0; 
      for (int j = 1; j <= p; ++j) {
        powers_ephi[j] = powers_ephi[j - 1] * ephi; 
      } 

      // Compute scaled modified spherical bessel function
      bessel_jn_scaled(p, omega * r, scale, bessel); 

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
    int p = builtin_helmholtz_table_->p(); 
    const double *sqf = builtin_helmholtz_table_->sqf(); 
    double omega = builtin_helmholtz_table_->omega(); 
    
    double *legendre = new double[(p + 1) * (p + 2) / 2]; 
    dcomplex_t *bessel = new dcomplex_t[p + 1]; 
    dcomplex_t *powers_ephi = new dcomplex_t[2 * p + 1]; 

    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, center); 
      double q = i->charge; 
      double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y()); 
      double r = dist.norm(); 

      // Compute cosine of the polar angle theta
      double ctheta = (r <= 1.0e-14 ? 1.0 : dist.z() / r); 

      // Compute exp(-i * phi) for the azimuthal angle phi 
      dcomplex_t ephi = (proj / r <= 1.0e-14 ? dcomplex_t{0.0, 0.0} : 
                         dcomplex_t{dist.x() / proj, -dist.y() / proj}); 

      // Compute powers of exp(-i * phi)
      powers_ephi[p] = 1.0; 
      for (int j = 1; j <= p; ++j) {
        powers_ephi[p + j] = powers_ephi[p + j - 1] * ephi; 
        powers_ephi[p - j] = conj(powers_ephi[p + j]);
      }
       
      // Compute scaled Hankel function 
      bessel_hn_scaled(p, omega * r, scale, bessel); 

      // Compute legendre polynomial 
      legendre_Plm(p, ctheta, legendre); 

      // Compute local expansion L_n^m 
      int curr = 0; 
      for (int n = 0; n <= p; ++n) {
        for (int m = -n; m <= n; ++m) {
          int idx = midx(n, fabs(m)); 
          L[curr++] += q * bessel[n] * sqf[idx] * legendre[idx] * 
            powers_ephi[p + m]; 
        }
      }
    }

    delete [] legendre; 
    delete [] bessel; 
    delete [] powers_ephi; 
    return std::unique_ptr<expansion_t>{retval}; 
  }

  std::unique_ptr<expansion_t> M_to_M(int from_child) const {
    /*
    // The function is called on the expansion of the child box and \p s_size is
    // the child box's size. 
    double h = s_size / 2; 
    Point center = views_.center(); 
    double px = center.x() + (from_child % 2 == 0 ? h : -h); 
    double py = center.y() + (from_child % 4 <= 1 ? h : -h); 
    double pz = center.z() + (from_child < 4 ? h : -h); 

    double scale = views_.scale(); 
    expansion_t *retval{new 
        expansion_t{Point{px, py, pz}, scale * 2, kSourcePrimary}}; 
    int p = builtin_helmholtz_table_->p(); 

    // Get precomputed Wigner d-matrix for rotation about the y-axis 
    const double *d1 = (from_child < 4 ? 
                        builtin_helmholtz_table_->dmat_plus(1.0 / sqrt(3.0)) : 
                        builtin_helmholtz_table_->dmat_plus(-1.0 / sqrt(3.0)));
    const double *d2 = (from_child < 4 ?
                        builtin_helmholtz_table_->dmat_minus(1.0 / sqrt(3.0)) : 
                        builtin_helmholtz_table_->dmat_minus(-1.0 / sqrt(3.0)));
    
    // Get precomputed coefficients for shifting along z-axis 
    const double *coeff = builtin_helmholtz_table_->m2m(scale); 

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

    rotate_sph_z(M, alpha, W1, false);  
    rotate_sph_y(W1, d1, W2, false); 

    for (int n = 0; n <= p; ++n) {
      for (int m = 0; m <= n; ++m) {
        dcomplex_t temp{0.0, 0.0}; 
        for (int np = m; np <= p; ++np) {
          temp += W2[midx(np, m)] * coeff[sidx(n, m, np, p)]; 
        }
        W1[midx(n, m)] = temp;
      }
    }

    rotate_sph_y(W1, d2, W2, false); 
    rotate_sph_z(W2, -alpha, W1, false); 

    delete [] W2; 
    return std::unique_ptr<expansion_t>{retval};
    */
    return std::unique_ptr<expansion_t>{nullptr};
  }

  std::unique_ptr<expansion_t> M_to_L(Index s_indx, Index t_index) const {
    return std::unique_ptr<expansion_t>{nullptr}; 
  }

  std::unique_ptr<expansion_t> L_to_L(int to_child, double t_size) const {
    // The function is called on the parent box and \p t_size is its child's
    // size. 
    Point center = views_.center(); 
    double h = t_size / 2; 
    double cx = center.x() + (to_child % 2 == 0 ? -h : h);
    double cy = center.y() + (to_child % 4 <= 1 ? -h : h);
    double cz = center.z() + (to_child < 4 ? -h : h);

    double scale = views_.scale();
    expansion_t *retval{new expansion_t{Point{cx, cy, cz},
          scale / 2, kTargetPrimary}};

    // Table of rotation angle about z-axis, as an integer multiple of pi / 4
    const int tab_alpha[8] = {1, 3, 7, 5, 1, 3, 7, 5};

    // Get rotation angle
    double alpha = tab_alpha[to_child] * M_PI_4;

    int p = builtin_helmholtz_table_->p();

    // Get precomputed Wigner d-matrix for rotation about the y-axis
    const double *d1 = (to_child < 4 ?
                        builtin_helmholtz_table_->dmat_plus(1.0 / sqrt(3)) :
                        builtin_helmholtz_table_->dmat_plus(-1.0 / sqrt(3)));
    const double *d2 = (to_child < 4 ?
                        builtin_helmholtz_table_->dmat_minus(1.0 / sqrt(3)) :
                        builtin_helmholtz_table_->dmat_minus(-1.0 / sqrt(3)));

    // Get precomputed coefficients for shifting along z-axis
    const double *coeff = builtin_helmholtz_table_->l2l(scale);

    // Get local expansion of the parent box
    dcomplex_t *L = reinterpret_cast<dcomplex_t *>(views_.view_data(0));

    // Temporary space for rotating local expansion
    dcomplex_t *W1 =
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    dcomplex_t *W2 = new dcomplex_t[(p + 1) * (p + 1)]; 

    rotate_sph_z(L, alpha, W1, true); 
    rotate_sph_y(W1, d1, W2, true); 

    for (int n = 0; n <= p; ++n) {
      for (int m = -n; m <= n; ++m) {
        dcomplex_t temp{0.0, 0.0}; 
        for (int np = fabs(m); np <= p; ++np) {
          temp += W2[lidx(np, m)] * coeff[sidx(n, fabs(m), np, p)]; 
        }
        W1[lidx(n, m)] = temp;
      }
    }

    rotate_sph_y(W1, d2, W2, true); 
    rotate_sph_z(W2, -alpha, W1, true); 

    delete [] W2; 
    return std::unique_ptr<expansion_t>{retval};
  }

  void M_to_T(Target *first, Target *last) const {
    int p = builtin_helmholtz_table_->p(); 
    double scale = views_.scale(); 
    double omega = builtin_helmholtz_table_->omega(); 
    double *legendre = new double[(p + 1) * (p + 2) / 2]; 
    dcomplex_t *bessel = new dcomplex_t[p + 1]; 
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
      dcomplex_t ephi = (proj / r <= 1.0e-14? dcomplex_t{1.0, 0.0} : 
                         dcomplex_t{dist.x() / proj, dist.y() / proj}); 

      // Compute powers of exp(i * phi)
      powers_ephi[0] = 1.0; 
      for (int j = 1; j <= p; ++j) {
        powers_ephi[j] = powers_ephi[j - 1] * ephi; 
      } 

      // Compute scaled hankel function 
      bessel_hn_scaled(p, omega * r, scale, bessel); 

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
    int p = builtin_helmholtz_table_->p(); 
    double scale = views_.scale(); 
    double omega = builtin_helmholtz_table_->omega(); 
    double *legendre = new double[(p + 1) * (p + 2) / 2]; 
    double *bessel = new double[p + 1]; 
    dcomplex_t *powers_ephi = new dcomplex_t[2 * p + 1]; 
    dcomplex_t *L = reinterpret_cast<dcomplex_t *>(views_.view_data(0)); 

    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, views_.center()); 
      dcomplex_t potential{0.0, 0.0}; 
      double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y()); 
      double r = dist.norm(); 

      // Compute cosine of the polar angle theta
      double ctheta = (r <= 1.0e-14 ? 1.0 : dist.z() / r); 
      
      // Compute exp(i * phi) for the azimuthal angle phi
      dcomplex_t ephi = (proj / r <= 1e-14? dcomplex_t{1.0, 0.0} : 
                         dcomplex_t{dist.x() / proj, dist.y() / proj}); 

      // Compute powers of exp(i * phi) 
      powers_ephi[p] = 1.0; 
      for (int j = 1; j <= p; ++j) {
        powers_ephi[p + j] = powers_ephi[p + j - 1] * ephi; 
        powers_ephi[p - j] = conj(powers_ephi[p + j]); 
      }

      // Compute scaled spherical Bessel function 
      bessel_jn_scaled(p, omega * r, scale, bessel); 

      // Compute legendre polynomial 
      legendre_Plm(p, ctheta, legendre); 

      // Evaluate L_n^m 
      int curr = 0; 
      for (int n = 0; n <= p; ++n) {
        for (int m = -n; m <= n; ++m) {
          potential += L[curr++] * powers_ephi[p + m] * bessel[n] * 
            legendre[midx(n, fabs(m))];
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
    double omega = builtin_helmholtz_table_->omega(); 
    for (auto i = t_first; i != t_last; ++i) {
      dcomplex_t potential{0.0, 0.0}; 
      for (auto j = s_first; j != s_last; ++j) {
        Point s2t = point_sub(i->position, j->position); 
        double dist = omega * s2t.norm(); 
        if (dist > 0) {
          dcomplex_t term{0.0, dist}; 
          potential += j->charge * exp(term) / term; 
        }
      }
      i->phi += potential;
    }
  }

  std::unique_ptr<expansion_t> M_to_I(Index s_index) const {
    double scale = views_.scale(); 
    expansion_t *retval{new expansion_t{views_.center(), 
          scale, kSourceIntermediate}}; 
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(views_.view_data(0)); 

    // Addresses of the views, in the order of x-, y-, and z-directions
    dcomplex_t *Prop[3] = {
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0)), 
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(3)),
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(6))}; 
    dcomplex_t *EvanP[3] = {
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(1)), 
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(4)),
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(7))}; 
    dcomplex_t *EvanM[3] = {
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(2)), 
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(4)),
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(8))}; 

    int p = builtin_helmholtz_table_->p(); 
    double wd = builtin_helmholtz_table_->omega() * 
      builtin_helmholtz_table_->size(scale); 

    // Get precomputed Wigner d-matrix 
    const double *d1 = builtin_helmholtz_table_->dmat_plus(0.0); 
    const double *d2 = builtin_helmholtz_table_->dmat_minus(0.0); 

    // Allocate temporary space to handle x-/y-direction expansion
    dcomplex_t *W1 = new dcomplex_t[(p + 1) * (p + 2) / 2]; 
    dcomplex_t *W2 = new dcomplex_t[(p + 1) * (p + 2) / 2]; 

    // Setup y-direction
    rotate_sph_z(M, -M_PI / 2, W1, false); 
    rotate_sph_y(W1, d2, W2, false); 

    // Setup x-direction
    rotate_sph_y(M, d1, W1, false); 

    // Addresses of the spherical harmonic expansions
    const dcomplex_t *SH[3] = {W1, W2, M}; 

    // Get quadratures 
    int s_e = builtin_helmholtz_table_->s_e(); 
    const double *x_e = builtin_helmholtz_table_->x_e(); 
    const int *f_e = builtin_helmholtz_table_->f_e(); 
    const int *m_e = builtin_helmholtz_table_->m_e(scale);
    const int *smf_e = builtin_helmholtz_table_->smf_e(scale); 
    const dcomplex_t *ealphaj_e = builtin_helmholtz_table_->ealphaj_e(scale); 

    int s_p = builtin_helmholtz_table_->s_p(); 
    const double *x_p = builtin_helmholtz_table_->x_p(); 
    const int *f_p = builtin_helmholtz_table_->f_p(); 
    const int *m_p = builtin_helmholtz_table_->m_p(scale); 
    const int *smf_p = builtin_helmholtz_table_->smf_p(scale); 
    const dcomplex_t *ealphaj_p = builtin_helmholtz_table_->ealphaj_p(scale); 

    double *legendre_e = new double[(p + 1) * (p + 2) / 2]; 
    dcomplex_t *legendre_p = new dcomplex_t[(p + 1) * (p + 2) / 2]; 

    for (int dir = 0; dir <= 2; ++dir) {
      int offset = 0; 

      // Compute propogating wave 
      for (int k = 0; k < s_p; ++k) {
        legendre_Plm_prop_scaled(p, cos(x_p[k]), scale, legendre_p); 

        dcomplex_t z0{0.0, 0.0}; 
        dcomplex_t *zp = new dcomplex_t[f_p[k] + 1](); // zp[0] not used
        dcomplex_t *zm = new dcomplex_t[f_p[k] + 1](); // zm[0] not used

        for (int n = 0; n <= p; ++n) {
          z0 += SH[dir][midx(n, 0)] * legendre_p[midx(n, 0)]; 
        }

        for (int m = 1; m <= f_p[k]; ++m) {
          for (int n = m; n <= p; ++n) {
            zp[m] += SH[dir][midx(n, m)] * legendre_p[midx(n, m)]; 
            zm[m] += conj(SH[dir][midx(n, m)]) * legendre_p[midx(n, m)]; 
          }
        }
        
        // Compute Prop(k, j) for 1 <= j <= m_p[k] / 2
        for (int j = 1; j <= m_p[k] / 2; ++j) {
          dcomplex_t temp = z0; 
          for (int m = 1; m <= f_p[k]; ++m) {
            int idx = smf_p[k] + (j - 1) * f_p[k] + m - 1; 
            temp += ealphaj_p[idx] * zp[m] + conj(ealphaj_p[idx]) * zm[m];
          }
          Prop[dir][offset++] = temp;
        }
        delete [] zp;
        delete [] zm; 
      }     

      // Compute evanescent wave
      offset = 0; 
      for (int k = 0; k < s_e; ++k) {
        legendre_Plm_evan_scaled(p, x_e[k] / wd, scale, legendre_e); 

        // Handle M_n^m where n is even 
        dcomplex_t *z1 = new dcomplex_t[f_e[k] + 1](); 
        // Handle M_n^m where n is odd 
        dcomplex_t *z2 = new dcomplex_t[f_e[k] + 1](); 

        for (int m = 0; m <= f_e[k]; m += 2) {
          for (int n = m; n <= p; n += 2) {
            z1[m] += SH[dir][midx(n, m)] * legendre_e[midx(n, m)]; 
          }

          for (int n = m + 1; n <= p; n += 2) {
            z2[m] += SH[dir][midx(n, m)] * legendre_e[midx(n, m)]; 
          }
        }

        for (int m = 1; m <= f_e[k]; m += 2) {
          for (int n = m; n <= p; n += 2) {
            z2[m] += SH[dir][midx(n, m)] * legendre_e[midx(n, m)]; 
          }
          for (int n = m + 1; n <= p; n += 2) {
            z1[m] += SH[dir][midx(n, m)] * legendre_e[midx(n, m)]; 
          }
        }

        // Compute Evan(k, j) for 1 <= j <= m_e[k] / 2
        for (int j = 1; j <= m_e[k] / 2; ++j) {
          dcomplex_t up{z1[0] + z2[0]}; 
          dcomplex_t dn{z1[0] - z2[0]}; 
          dcomplex_t power_I{0.0, 1.0}; 
          for (int  m = 1; m <= f_e[k]; ++m) {
            int idx = smf_e[k] + (j - 1) * f_e[k] + m - 1;  
            up += 2 * real(ealphaj_e[idx] * (z1[m] + z2[m])) * power_I; 
            dn += 2 * real(ealphaj_e[idx] * (z1[m] - z2[m])) * power_I; 
            power_I *= dcomplex_t{0.0, 1.0}; 
          }
          EvanP[dir][offset] = up; 
          EvanM[dir][offset] = dn; 
          offset++; 
        }
        
        delete [] z1; 
        delete [] z2; 
      }
    }

    delete [] W1; 
    delete [] W2; 
    delete [] legendre_p; 
    delete [] legendre_e;  
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
    
    // Addresses of the views on the source side
    dcomplex_t *Prop_x = 
      reinterpret_cast<dcomplex_t *>(views_.view_data(0)); 
    dcomplex_t *Evan_px = 
      reinterpret_cast<dcomplex_t *>(views_.view_data(1)); 
    dcomplex_t *Evan_mx = 
      reinterpret_cast<dcomplex_t *>(views_.view_data(2)); 
    dcomplex_t *Prop_y = 
      reinterpret_cast<dcomplex_t *>(views_.view_data(3)); 
    dcomplex_t *Evan_py = 
      reinterpret_cast<dcomplex_t *>(views_.view_data(4)); 
    dcomplex_t *Evan_my = 
      reinterpret_cast<dcomplex_t *>(views_.view_data(5)); 
    dcomplex_t *Prop_z = 
      reinterpret_cast<dcomplex_t *>(views_.view_data(6)); 
    dcomplex_t *Evan_pz = 
      reinterpret_cast<dcomplex_t *>(views_.view_data(7)); 
    dcomplex_t *Evan_mz = 
      reinterpret_cast<dcomplex_t *>(views_.view_data(8)); 

    ViewSet views{kTargetIntermediate, Point{px, py, pz}, 2 * scale}; 

    // Each S will generate from 1 to 3 views (evan + prop) on the target
    // side. For propagating wave, the terms are doubled as conjugacy is lost
    // after shifting
    int n_e = builtin_helmholtz_table_->n_e(scale); 
    int n_p = builtin_helmholtz_table_->n_p(scale) * 2; 
    size_t bytes_e = n_e * sizeof(dcomplex_t); 
    size_t bytes_p = n_p * sizeof(dcomplex_t); 
    
    dcomplex_t *T1 = new dcomplex_t[n_p](); 
    dcomplex_t *T2 = new dcomplex_t[n_e](); 
    dcomplex_t *T3 = new dcomplex_t[n_p](); 
    dcomplex_t *T4 = new dcomplex_t[n_e](); 
    dcomplex_t *T5 = new dcomplex_t[n_p](); 
    dcomplex_t *T6 = new dcomplex_t[n_e](); 
    char *C1 = reinterpret_cast<char *>(T1);
    char *C2 = reinterpret_cast<char *>(T2);
    char *C3 = reinterpret_cast<char *>(T3);
    char *C4 = reinterpret_cast<char *>(T4);
    char *C5 = reinterpret_cast<char *>(T5);
    char *C6 = reinterpret_cast<char *>(T6);
    dcomplex_t *T[6] = {T1, T2, T3, T4, T5, T6}; 
    char *C[6] = {C1, C2, C3, C4, C5, C6}; 
    bool used[3] = {false, false, false}; 

    for (int i = 0; i < 3; ++i) {
      int j = 2 * i; 
      int tag = merge_and_shift_table[dx + 2][dy + 2][dz + 2][i];

      if (tag == -1) {
        break;
      }

      if (tag <= 1) {
        e2e_p(T[j], Prop_z, dx, dy, 0, scale, true, false); 
        e2e_e(T[j + 1], Evan_mz, dx, dy, 0, scale); 
      } else if (tag <= 5) {
        e2e_p(T[j], Prop_y, dx, dy, 0, scale, true, false); 
        e2e_e(T[j + 1], Evan_my, dz, dx, 0, scale);
      } else if (tag <= 13) {
        e2e_p(T[j], Prop_x, -dz, dy, 0, scale, true, false); 
        e2e_e(T[j + 1], Evan_mx, -dz, dy, 0, scale); 
      } else if (tag <= 15) {
        e2e_p(T[j], Prop_z, -dx, -dy, 0, scale, false, false);
        e2e_e(T[j + 1], Evan_pz, -dx, -dy, 0, scale);
      } else if (tag <= 19) {
        e2e_p(T[j], Prop_y, -dz, -dx, 0, scale, false, false);
        e2e_e(T[j + 1], Evan_py, -dz, -dx, 0, scale);
      } else {
        e2e_p(T[j], Prop_x, dz, -dy, 0, scale, false, false); 
        e2e_e(T[j + 1], Evan_px, dz, -dy, 0, scale);
      }

      views.add_view(2 * tag, bytes_p, C[j]); 
      views.add_view(2 * tag + 1, bytes_e, C[j + 1]); 
      used[i] = true;
    }

    if (used[1] == false) {
      delete [] T3;
      delete [] T4; 
    }

    if (used[2] == false) {
      delete [] T5;
      delete [] T6; 
    }

    expansion_t *retval = new expansion_t{views};
    return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> I_to_L(Index t_index, double t_size) const {
    // t_index and t_size are the index_and size of the child
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
    dcomplex_t *Evan[28]{nullptr}; 
    dcomplex_t *Prop[28]{nullptr}; 
    for (int i = 0; i < 28; ++i) {
      Prop[i] = reinterpret_cast<dcomplex_t *>(views_.view_data(2 * i));
      Evan[i] = reinterpret_cast<dcomplex_t *>(views_.view_data(2 * i + 1)); 
    }
    dcomplex_t *L = 
      reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    
    int n_e = builtin_helmholtz_table_->n_e(scale); 
    int n_p = builtin_helmholtz_table_->n_p(scale) * 2; 
    dcomplex_t *S = new dcomplex_t[(n_e + n_p) * 6](); 
    dcomplex_t *sProp_mz = S; 
    dcomplex_t *sProp_pz = sProp_mz + n_p; 
    dcomplex_t *sProp_my = sProp_pz + n_p; 
    dcomplex_t *sProp_py = sProp_my + n_p; 
    dcomplex_t *sProp_mx = sProp_py + n_p;
    dcomplex_t *sProp_px = sProp_mx + n_p; 
    dcomplex_t *sEvan_mz = sProp_px + n_p;
    dcomplex_t *sEvan_pz = sEvan_mz + n_e; 
    dcomplex_t *sEvan_my = sEvan_pz + n_e;
    dcomplex_t *sEvan_py = sEvan_my + n_e; 
    dcomplex_t *sEvan_mx = sEvan_py + n_e; 
    dcomplex_t *sEvan_px = sEvan_mx + n_e; 

    switch (to_child) {
    case 0:
      e2e_e(sEvan_mz, Evan[uall], 0, 0, 3, scale);
      e2e_e(sEvan_mz, Evan[u1234], 0, 0, 2, scale);
      e2e_e(sEvan_pz, Evan[dall], 0, 0, 2, scale);

      e2e_e(sEvan_my, Evan[nall], 0, 0, 3, scale);
      e2e_e(sEvan_my, Evan[n1256], 0, 0, 2, scale);
      e2e_e(sEvan_my, Evan[n12], 0, 0, 2, scale);
      e2e_e(sEvan_py, Evan[sall], 0, 0, 2, scale);

      e2e_e(sEvan_mx, Evan[eall], 0, 0, 3, scale);
      e2e_e(sEvan_mx, Evan[e1357], 0, 0, 2, scale);
      e2e_e(sEvan_mx, Evan[e13], 0, 0, 2, scale);
      e2e_e(sEvan_mx, Evan[e1], 0, 0, 2, scale);
      e2e_e(sEvan_px, Evan[wall], 0, 0, 2, scale);

      e2e_p(sProp_mz, Prop[uall], 0, 0, 3, scale, false, true);
      e2e_p(sProp_mz, Prop[u1234], 0, 0, 2, scale, false, true);
      e2e_p(sProp_pz, Prop[dall], 0, 0, 2, scale, false, true);

      e2e_p(sProp_my, Prop[nall], 0, 0, 3, scale, false, true);
      e2e_p(sProp_my, Prop[n1256], 0, 0, 2, scale, false, true);
      e2e_p(sProp_my, Prop[n12], 0, 0, 2, scale, false, true);
      e2e_p(sProp_py, Prop[sall], 0, 0, 2, scale, false, true);

      e2e_p(sProp_mx, Prop[eall], 0, 0, 3, scale, false, true);
      e2e_p(sProp_mx, Prop[e1357], 0, 0, 2, scale, false, true);
      e2e_p(sProp_mx, Prop[e13], 0, 0, 2, scale, false, true);
      e2e_p(sProp_mx, Prop[e1], 0, 0, 2, scale, false, true);
      e2e_p(sProp_px, Prop[wall], 0, 0, 2, scale, false, true);
      break;
    case 1:
      e2e_e(sEvan_mz, Evan[uall], -1, 0, 3, scale);
      e2e_e(sEvan_mz, Evan[u1234], -1, 0, 2, scale);
      e2e_e(sEvan_pz, Evan[dall], 1, 0, 2, scale);

      e2e_e(sEvan_my, Evan[nall], 0, -1, 3, scale);
      e2e_e(sEvan_my, Evan[n1256], 0, -1, 2, scale);
      e2e_e(sEvan_my, Evan[n12], 0, -1, 2, scale);
      e2e_e(sEvan_py, Evan[sall], 0, 1, 2, scale);

      e2e_e(sEvan_mx, Evan[eall], 0, 0, 2, scale);
      e2e_e(sEvan_px, Evan[wall], 0, 0, 3, scale);
      e2e_e(sEvan_px, Evan[w2468], 0, 0, 2, scale);
      e2e_e(sEvan_px, Evan[w24], 0, 0, 2, scale);
      e2e_e(sEvan_px, Evan[w2], 0, 0, 2, scale);

      e2e_p(sProp_mz, Prop[uall], -1, 0, 3, scale, false, true);
      e2e_p(sProp_mz, Prop[u1234], -1, 0, 2, scale, false, true);
      e2e_p(sProp_pz, Prop[dall], 1, 0, 2, scale, false, true);

      e2e_p(sProp_my, Prop[nall], 0, -1, 3, scale, false, true);
      e2e_p(sProp_my, Prop[n1256], 0, -1, 2, scale, false, true);
      e2e_p(sProp_my, Prop[n12], 0, -1, 2, scale, false, true);
      e2e_p(sProp_py, Prop[sall], 0, 1, 2, scale, false, true);

      e2e_p(sProp_mx, Prop[eall], 0, 0, 2, scale, false, true);
      e2e_p(sProp_px, Prop[wall], 0, 0, 3, scale, false, true);
      e2e_p(sProp_px, Prop[w2468], 0, 0, 2, scale, false, true);
      e2e_p(sProp_px, Prop[w24], 0, 0, 2, scale, false, true);
      e2e_p(sProp_px, Prop[w2], 0, 0, 2, scale, false, true);
      break;
    case 2:
      e2e_e(sEvan_mz, Evan[uall], 0, -1, 3, scale);
      e2e_e(sEvan_mz, Evan[u1234], 0, -1, 2, scale);
      e2e_e(sEvan_pz, Evan[dall], 0, 1, 2, scale);

      e2e_e(sEvan_my, Evan[nall], 0, 0, 2, scale);
      e2e_e(sEvan_py, Evan[sall], 0, 0, 3, scale);
      e2e_e(sEvan_py, Evan[s3478], 0, 0, 2, scale);
      e2e_e(sEvan_py, Evan[s34], 0, 0, 2, scale);

      e2e_e(sEvan_mx, Evan[eall], 0, -1, 3, scale);
      e2e_e(sEvan_mx, Evan[e1357], 0, -1, 2, scale);
      e2e_e(sEvan_mx, Evan[e13], 0, -1, 2, scale);
      e2e_e(sEvan_mx, Evan[e3], 0, -1, 2, scale);
      e2e_e(sEvan_px, Evan[wall], 0, 1, 2, scale);

      e2e_p(sProp_mz, Prop[uall], 0, -1, 3, scale, false, true);
      e2e_p(sProp_mz, Prop[u1234], 0, -1, 2, scale, false, true);
      e2e_p(sProp_pz, Prop[dall], 0, 1, 2, scale, false, true);

      e2e_p(sProp_my, Prop[nall], 0, 0, 2, scale, false, true);
      e2e_p(sProp_py, Prop[sall], 0, 0, 3, scale, false, true);
      e2e_p(sProp_py, Prop[s3478], 0, 0, 2, scale, false, true);
      e2e_p(sProp_py, Prop[s34], 0, 0, 2, scale, false, true);

      e2e_p(sProp_mx, Prop[eall], 0, -1, 3, scale, false, true);
      e2e_p(sProp_mx, Prop[e1357], 0, -1, 2, scale, false, true);
      e2e_p(sProp_mx, Prop[e13], 0, -1, 2, scale, false, true);
      e2e_p(sProp_mx, Prop[e3], 0, -1, 2, scale, false, true);
      e2e_p(sProp_px, Prop[wall], 0, 1, 2, scale, false, true);
      break;
    case 3:
      e2e_e(sEvan_mz, Evan[uall], -1, -1, 3, scale);
      e2e_e(sEvan_mz, Evan[u1234], -1, -1, 2, scale);
      e2e_e(sEvan_pz, Evan[dall], 1, 1, 2, scale);

      e2e_e(sEvan_my, Evan[nall], 0, -1, 2, scale);
      e2e_e(sEvan_py, Evan[sall], 0, 1, 3, scale);
      e2e_e(sEvan_py, Evan[s3478], 0, 1, 2, scale);
      e2e_e(sEvan_py, Evan[s34], 0, 1, 2, scale);

      e2e_e(sEvan_mx, Evan[eall], 0, -1, 2, scale);
      e2e_e(sEvan_px, Evan[wall], 0, 1, 3, scale);
      e2e_e(sEvan_px, Evan[w2468], 0, 1, 2, scale);
      e2e_e(sEvan_px, Evan[w24], 0, 1, 2, scale);
      e2e_e(sEvan_px, Evan[w4], 0, 1, 2, scale);

      e2e_p(sProp_mz, Prop[uall], -1, -1, 3, scale, false, true);
      e2e_p(sProp_mz, Prop[u1234], -1, -1, 2, scale, false, true);
      e2e_p(sProp_pz, Prop[dall], 1, 1, 2, scale, false, true);

      e2e_p(sProp_my, Prop[nall], 0, -1, 2, scale, false, true);
      e2e_p(sProp_py, Prop[sall], 0, 1, 3, scale, false, true);
      e2e_p(sProp_py, Prop[s3478], 0, 1, 2, scale, false, true);
      e2e_p(sProp_py, Prop[s34], 0, 1, 2, scale, false, true);

      e2e_p(sProp_mx, Prop[eall], 0, -1, 2, scale, false, true);
      e2e_p(sProp_px, Prop[wall], 0, 1, 3, scale, false, true);
      e2e_p(sProp_px, Prop[w2468], 0, 1, 2, scale, false, true);
      e2e_p(sProp_px, Prop[w24], 0, 1, 2, scale, false, true);
      e2e_p(sProp_px, Prop[w4], 0, 1, 2, scale, false, true);
      break;
    case 4:
      e2e_e(sEvan_mz, Evan[uall], 0, 0, 2, scale);
      e2e_e(sEvan_pz, Evan[dall], 0, 0, 3, scale);
      e2e_e(sEvan_pz, Evan[d5678], 0, 0, 2, scale);

      e2e_e(sEvan_my, Evan[nall], -1, 0, 3, scale);
      e2e_e(sEvan_my, Evan[n1256], -1, 0, 2, scale);
      e2e_e(sEvan_my, Evan[n56], -1, 0, 2, scale);
      e2e_e(sEvan_py, Evan[sall], 1, 0, 2, scale);

      e2e_e(sEvan_mx, Evan[eall], 1, 0, 3, scale);
      e2e_e(sEvan_mx, Evan[e1357], 1, 0, 2, scale);
      e2e_e(sEvan_mx, Evan[e57], 1, 0, 2, scale);
      e2e_e(sEvan_mx, Evan[e5], 1, 0, 2, scale);
      e2e_e(sEvan_px, Evan[wall], -1, 0, 2, scale);

      e2e_p(sProp_mz, Prop[uall], 0, 0, 2, scale, false, true);
      e2e_p(sProp_pz, Prop[dall], 0, 0, 3, scale, false, true);
      e2e_p(sProp_pz, Prop[d5678], 0, 0, 2, scale, false, true);

      e2e_p(sProp_my, Prop[nall], -1, 0, 3, scale, false, true);
      e2e_p(sProp_my, Prop[n1256], -1, 0, 2, scale, false, true);
      e2e_p(sProp_my, Prop[n56], -1, 0, 2, scale, false, true);
      e2e_p(sProp_py, Prop[sall], 1, 0, 2, scale, false, true);

      e2e_p(sProp_mx, Prop[eall], 1, 0, 3, scale, false, true);
      e2e_p(sProp_mx, Prop[e1357], 1, 0, 2, scale, false, true);
      e2e_p(sProp_mx, Prop[e57], 1, 0, 2, scale, false, true);
      e2e_p(sProp_mx, Prop[e5], 1, 0, 2, scale, false, true);
      e2e_p(sProp_px, Prop[wall], -1, 0, 2, scale, false, true);
      break;
    case 5:
      e2e_e(sEvan_mz, Evan[uall], -1, 0, 2, scale);
      e2e_e(sEvan_pz, Evan[dall], 1, 0, 3, scale);
      e2e_e(sEvan_pz, Evan[d5678], 1, 0, 2, scale);

      e2e_e(sEvan_my, Evan[nall], -1, -1, 3, scale);
      e2e_e(sEvan_my, Evan[n1256], -1, -1, 2, scale);
      e2e_e(sEvan_my, Evan[n56], -1, -1, 2, scale);
      e2e_e(sEvan_py, Evan[sall], 1, 1, 2, scale);

      e2e_e(sEvan_mx, Evan[eall], 1, 0, 2, scale);
      e2e_e(sEvan_px, Evan[wall], -1, 0, 3, scale);
      e2e_e(sEvan_px, Evan[w2468], -1, 0, 2, scale);
      e2e_e(sEvan_px, Evan[w68], -1, 0, 2, scale);
      e2e_e(sEvan_px, Evan[w6], -1, 0, 2, scale);

      e2e_p(sProp_mz, Prop[uall], -1, 0, 2, scale, false, true);
      e2e_p(sProp_pz, Prop[dall], 1, 0, 3, scale, false, true);
      e2e_p(sProp_pz, Prop[d5678], 1, 0, 2, scale, false, true);

      e2e_p(sProp_my, Prop[nall], -1, -1, 3, scale, false, true);
      e2e_p(sProp_my, Prop[n1256], -1, -1, 2, scale, false, true);
      e2e_p(sProp_my, Prop[n56], -1, -1, 2, scale, false, true);
      e2e_p(sProp_py, Prop[sall], 1, 1, 2, scale, false, true);

      e2e_p(sProp_mx, Prop[eall], 1, 0, 2, scale, false, true);
      e2e_p(sProp_px, Prop[wall], -1, 0, 3, scale, false, true);
      e2e_p(sProp_px, Prop[w2468], -1, 0, 2, scale, false, true);
      e2e_p(sProp_px, Prop[w68], -1, 0, 2, scale, false, true);
      e2e_p(sProp_px, Prop[w6], -1, 0, 2, scale, false, true);
      break;
    case 6:
      e2e_e(sEvan_mz, Evan[uall], 0, -1, 2, scale);
      e2e_e(sEvan_pz, Evan[dall], 0, 1, 3, scale);
      e2e_e(sEvan_pz, Evan[d5678], 0, 1, 2, scale);

      e2e_e(sEvan_my, Evan[nall], -1, 0, 2, scale);
      e2e_e(sEvan_py, Evan[sall], 1, 0, 3, scale);
      e2e_e(sEvan_py, Evan[s3478], 1, 0, 2, scale);
      e2e_e(sEvan_py, Evan[s78], 1, 0, 2, scale);

      e2e_e(sEvan_mx, Evan[eall], 1, -1, 3, scale);
      e2e_e(sEvan_mx, Evan[e1357], 1, -1, 2, scale);
      e2e_e(sEvan_mx, Evan[e57], 1, -1, 2, scale);
      e2e_e(sEvan_mx, Evan[e7], 1, -1, 2, scale);
      e2e_e(sEvan_px, Evan[wall], -1, 1, 2, scale);

      e2e_p(sProp_mz, Prop[uall], 0, -1, 2, scale, false, true);
      e2e_p(sProp_pz, Prop[dall], 0, 1, 3, scale, false, true);
      e2e_p(sProp_pz, Prop[d5678], 0, 1, 2, scale, false, true);

      e2e_p(sProp_my, Prop[nall], -1, 0, 2, scale, false, true);
      e2e_p(sProp_py, Prop[sall], 1, 0, 3, scale, false, true);
      e2e_p(sProp_py, Prop[s3478], 1, 0, 2, scale, false, true);
      e2e_p(sProp_py, Prop[s78], 1, 0, 2, scale, false, true);

      e2e_p(sProp_mx, Prop[eall], 1, -1, 3, scale, false, true);
      e2e_p(sProp_mx, Prop[e1357], 1, -1, 2, scale, false, true);
      e2e_p(sProp_mx, Prop[e57], 1, -1, 2, scale, false, true);
      e2e_p(sProp_mx, Prop[e7], 1, -1, 2, scale, false, true);
      e2e_p(sProp_px, Prop[wall], -1, 1, 2, scale, false, true);
      break;
    case 7:
      e2e_e(sEvan_mz, Evan[uall], -1, -1, 2, scale);
      e2e_e(sEvan_pz, Evan[dall], 1, 1, 3, scale);
      e2e_e(sEvan_pz, Evan[d5678], 1, 1, 2, scale);

      e2e_e(sEvan_my, Evan[nall], -1, -1, 2, scale);
      e2e_e(sEvan_py, Evan[sall], 1, 1, 3, scale);
      e2e_e(sEvan_py, Evan[s3478], 1, 1, 2, scale);
      e2e_e(sEvan_py, Evan[s78], 1, 1, 2, scale);

      e2e_e(sEvan_mx, Evan[eall], 1, -1, 2, scale);
      e2e_e(sEvan_px, Evan[wall], -1, 1, 3, scale);
      e2e_e(sEvan_px, Evan[w2468], -1, 1, 2, scale);
      e2e_e(sEvan_px, Evan[w68], -1, 1, 2, scale);
      e2e_e(sEvan_px, Evan[w8], -1, 1, 2, scale);

      e2e_p(sProp_mz, Prop[uall], -1, -1, 2, scale, false, true);
      e2e_p(sProp_pz, Prop[dall], 1, 1, 3, scale, false, true);
      e2e_p(sProp_pz, Prop[d5678], 1, 1, 2, scale, false, true);

      e2e_p(sProp_my, Prop[nall], -1, -1, 2, scale, false, true);
      e2e_p(sProp_py, Prop[sall], 1, 1, 3, scale, false, true);
      e2e_p(sProp_py, Prop[s3478], 1, 1, 2, scale, false, true);
      e2e_p(sProp_py, Prop[s78], 1, 1, 2, scale, false, true);

      e2e_p(sProp_mx, Prop[eall], 1, -1, 2, scale, false, true);
      e2e_p(sProp_px, Prop[wall], -1, 1, 3, scale, false, true);
      e2e_p(sProp_px, Prop[w2468], -1, 1, 2, scale, false, true);
      e2e_p(sProp_px, Prop[w68], -1, 1, 2, scale, false, true);
      e2e_p(sProp_px, Prop[w8], -1, 1, 2, scale, false, true);
      break;
    }

    e2l(sProp_mz, sEvan_mz, 'z', false, L); 
    e2l(sProp_pz, sEvan_pz, 'z', true, L); 
    e2l(sProp_my, sEvan_my, 'y', false, L); 
    e2l(sProp_py, sEvan_py, 'y', true, L);
    e2l(sProp_mx, sEvan_mx, 'x', false, L); 
    e2l(sProp_px, sEvan_px, 'x', true, L); 

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
    update_helmholtz_table(n_digits, domain_size, kernel_params[0]); 
  }

  static void delete_table() { }

  static double compute_scale(Index index) {
    return builtin_helmholtz_table_->scale(index.level()); 
  }

  static int weight_estimate(Operation op,
                             Index s = Index{}, Index t = Index{}) {
    // The weight needs to updated 
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

  void rotate_sph_z(const dcomplex_t *M, double alpha, dcomplex_t *MR, 
                    bool is_local) const {
    int p = builtin_helmholtz_table_->p(); 
    // Compute exp(i * alpha)
    dcomplex_t ealpha = dcomplex_t{cos(alpha), sin(alpha)}; 

    // Compute powers of exp(i * alpha)
    dcomplex_t *powers_ealpha = new dcomplex_t[2 * p + 1]; 
    powers_ealpha[p] = dcomplex_t{1.0, 0.0}; 
    for (int j = 1; j <= p; ++j) {
      powers_ealpha[p + j] = powers_ealpha[p + j - 1] * ealpha; 
      powers_ealpha[p - j] = conj(powers_ealpha[p + j]); 
    }

    int offset = 0; 
    for (int n = 0; n <= p; ++n) {
      int mmin = (is_local ? -n : 0); 
      for (int m = mmin; m <= n; ++m) {
        MR[offset] = M[offset] * powers_ealpha[p + m]; 
        offset++;
      }
    }
  }

  void rotate_sph_y(const dcomplex_t *M, const double *d, dcomplex_t *MR, 
                    bool is_local) const {
    int p = builtin_helmholtz_table_->p(); 
    int curr = 0; 

    if (is_local) {
      for (int n = 0; n <= p; ++n) {
        for (int mp = -n; mp <= -1; mp++) {
          // Get L_n^0 
          const dcomplex_t *Ln = &M[lidx(n, 0)]; 
          // Apply d(n, mp, m) = d(n, -m, -mp) 
          for (int m = -n; m <= -1; ++m) 
            MR[curr] += Ln[m] * d[didx(n, -m, -mp)]; 
          // Apply d(n, mp, m) = (-1)^(m - mp) d(n, m, mp)
          // and the input M needs to be scaled (-1)^m 
          for (int m = 0; m <= n; ++m) 
            MR[curr] += Ln[m] * d[didx(n, m, mp)] * pow(-1, mp);
          curr++; 
        } 

        int power_mp = 1; 
        for (int mp = 0; mp <= n; ++mp) {
          // Get d_n^{mp, 0}
          const double *coeff = &d[didx(n, mp, 0)]; 
          // Get L_n^0 
          const dcomplex_t *Ln = &M[lidx(n, 0)]; 
          MR[curr] = Ln[0] * coeff[0]; 
          double power_m = -1; 
          for (int m = 1; m <= n; ++m) {
            MR[curr] += (Ln[m] * power_m * coeff[m] + Ln[-m] * coeff[-m]);
            power_m = -power_m;
          }
          MR[curr++] *= power_mp; 
          power_mp = -power_mp;
        }
      }
    } else {
      for (int n = 0; n <= p; ++n) {
        int power_mp = 1; 
        for (int mp = 0; mp <= n; ++mp) {
          // Get d_n^{mp, 0}
          const double *coeff = &d[didx(n, mp, 0)]; 
          // Get M_n^0 
          const dcomplex_t *Mn = &M[midx(n, 0)]; 
          MR[curr] = Mn[0] * coeff[0]; 
          double power_m = -1; 
          for (int m = 1; m <= n; ++m) {
            MR[curr] += (Mn[m] * power_m * coeff[m] + conj(Mn[m]) * coeff[-m]); 
            power_m = -power_m;
          }
          MR[curr++] *= power_mp; 
          power_mp = -power_mp;
        }
      }
    }
  }

  void e2e_p(dcomplex_t *M, const dcomplex_t *W, int x, int y, int z, 
             double scale, bool conjugate, bool target) const {
    const dcomplex_t *xs = builtin_helmholtz_table_->xs_p(scale); 
    const dcomplex_t *ys = builtin_helmholtz_table_->ys_p(scale); 
    const dcomplex_t *zs = builtin_helmholtz_table_->zs_p(scale); 
    const int *m = builtin_helmholtz_table_->m_e(scale); 
    const int *sm = builtin_helmholtz_table_->sm_e(scale); 
    int s = builtin_helmholtz_table_->s_e(); 

    // There are two cases handled by this function. 
    // (1) Shifting propagating wave from source side to target side
    // (2) Shifting propagating wave from target side to target side 
    // For case 1, the input is saved for W(k, j), 1 <= j <= m[k] / 2
    // and the output is saved for 1 <= j <= m[k]
    // For case 2, both input and output are saved for 1 <= j <= m[k]

    if (target) {
      for (int k = 0; k < s; ++k) {
        int curr = sm[k]; 
        dcomplex_t factor_z = zs[7 * k + 3 + z]; 
        int mk2 = m[k] / 2; 
        for (int j = 0; j < m[k] / 2; ++j) {
          int offset = (curr + j) * 7 + 3; 
          dcomplex_t factor_xy = xs[offset + x] * ys[offset + y]; 
          M[2 * curr + j] += W[2 * curr + j] * factor_z * factor_xy; 
          M[2 * curr + j + mk2] += 
            W[2 * curr + j + mk2] * factor_z * conj(factor_xy); 
        }
      }
    } else {
      if (conjugate) {
        for (int k = 0; k < s; ++k) {
          int curr = sm[k]; 
          dcomplex_t factor_z = zs[7 * k + 3 + z]; 
          int mk2 = m[k] / 2; 
          for (int j = 0; j < m[k] / 2; ++j) {
            int offset = (curr + j) * 7 + 3; 
            dcomplex_t factor_xy = xs[offset + x] * ys[offset+y]; 
            M[2 * curr + j] += conj(W[curr + j]) * factor_z * factor_xy; 
            M[2 * curr + j + mk2] += W[curr + j] * factor_z * conj(factor_xy); 
          }
        }
      } else {
        for (int k = 0; k < s; ++k) {
          int curr = sm[k]; 
          dcomplex_t factor_z = zs[7 * k + 3 + z]; 
          int mk2 = m[k] / 2; 
          for (int j = 0; j < m[k] / 2; ++j) {
            int offset = (curr + j) * 7 + 3; 
            dcomplex_t factor_xy = xs[offset + x] * ys[offset + y]; 
            M[2 * curr + j] += W[curr + j] * factor_z * factor_xy;
            M[2 * curr + j + mk2] += conj(W[curr + j] * factor_xy) * factor_z; 
          }
        }
      }
    }
  }
  
  void e2e_e(dcomplex_t *M, const dcomplex_t *W, int x, int y, int z, 
             double scale) const {
    const dcomplex_t *xs = builtin_helmholtz_table_->xs_e(scale); 
    const dcomplex_t *ys = builtin_helmholtz_table_->ys_e(scale); 
    const double *zs = builtin_helmholtz_table_->zs_e(scale); 
    const int *m = builtin_helmholtz_table_->m_e(scale); 
    const int *sm = builtin_helmholtz_table_->sm_e(scale); 
    int s = builtin_helmholtz_table_->s_e(); 

    for (int k = 0; k < s; ++k) {
      double factor_z = zs[4 * k + z]; 
      for (int j = 0; j < m[k] / 2; ++j) {
        int offset = (sm[k] + j) * 7 + 3; 
        dcomplex_t factor_xy = xs[offset + x] * ys[offset + y]; 
        M[offset] += W[offset] * factor_z * factor_xy; 
        offset++;
      }
    }
  }
  
  void e2l(const dcomplex_t *Prop, const dcomplex_t *Evan, char dir, 
           bool sgn, dcomplex_t *L) const {
    // Note: this function is called on the parent node 
    double scale = views_.scale() / 2; 
    int p = builtin_helmholtz_table_->p(); 
    const double *sqf = builtin_helmholtz_table_->sqf(); 

    double *legendre_e = new double[(p + 1) * (p + 2) / 2]; 
    int s_e = builtin_helmholtz_table_->s_e(); 
    const int *m_e = builtin_helmholtz_table_->m_e(scale); 
    const int *sm_e = builtin_helmholtz_table_->sm_e(scale); 
    const int *smf_e = builtin_helmholtz_table_->smf_e(scale); 
    const int *f_e = builtin_helmholtz_table_->f_e(); 
    const dcomplex_t *ealphaj_e = builtin_helmholtz_table_->ealphaj_e(scale); 
    const double *x_e = builtin_helmholtz_table_->x_e(); 
    const double *w_e = builtin_helmholtz_table_->w_e(); 
    double wd = builtin_helmholtz_table_->omega() * 
      builtin_helmholtz_table_->size(scale); 

    dcomplex_t *legendre_p = new dcomplex_t[(p + 1) * (p + 2) / 2]; 
    int s_p = builtin_helmholtz_table_->s_p();
    const int *m_p = builtin_helmholtz_table_->m_p(scale); 
    const int *sm_p = builtin_helmholtz_table_->sm_p(scale); 
    const int *smf_p = builtin_helmholtz_table_->smf_p(scale); 
    const int *f_p = builtin_helmholtz_table_->f_p(); 
    const dcomplex_t *ealphaj_p = builtin_helmholtz_table_->ealphaj_p(scale); 
    const double *x_p = builtin_helmholtz_table_->x_p(); 
    const double *w_p = builtin_helmholtz_table_->w_p(); 

    dcomplex_t *contrib = nullptr; 
    dcomplex_t *W1 = new dcomplex_t[(p + 1) * (p + 1)]; 
    dcomplex_t *W2 = new dcomplex_t[(p + 1) * (p + 1)]; 

    // Convert evanescent wave into local expansion
    for (int k = 0; k < s_e; ++k) {
      int mk2 = m_e[k] / 2; 

      // Note: Evan(k, j) needs to be scaled by -i * w_evan(k) / wd / m_evan(k). 
      // Without this factor and if m_evan[k] is even
      // Evan(k, j) and Evan(k, j + mk2) are conjugate of each other      
      
      // Computes sum_{j=1}^{m_evan(k)} Evan(k, j) e^{-i * m * alpha_j}      
      dcomplex_t *z = new dcomplex_t[f_e[k] + 1](); 

      // m = 0
      for (int j = 1; j <= mk2; ++j) {
        int idx = sm_e[k] + j - 1; 
        z[0] += 2 * real(Evan[idx]); 
      }

      // m = 1, ..., f_evan[k], where m is odd 
      for (int m = 1; m <= f_e[k]; m += 2) {
        for (int j = 1; j <= mk2; ++j) {
          int widx = sm_e[k] + j - 1; 
          int aidx = smf_e[k] + (j - 1) * f_e[k] + m - 1; 
          z[m] += conj(ealphaj_e[aidx]) * 2.0 * imag(Evan[widx]); 
        }
      }

      // m = 2, ..., f_evan[k], where m is even 
      for (int m = 2; m <= f_e[k]; m += 2) {
        for (int j = 1; j <= mk2; ++j) {
          int widx = sm_e[k] + j - 1; 
          int aidx = smf_e[k] + (j - 1) * f_e[k] + m - 1; 
          z[m] += conj(ealphaj_e[aidx]) * 2.0 * real(Evan[widx]);
        }
      }

      legendre_Plm_evan_scaled(p, x_e[k] / wd, scale, legendre_e); 
      
      double factor1 = 2 * w_e[k] / wd / m_e[k]; 
      for (int n = 0; n <= p; ++n) {
        int mmax = (n <= f_e[k] ? n : f_e[k]); 
        
        // Process L_n^0
        W1[lidx(n, 0)] += z[0] * legendre_e[midx(n, 0)] * 
          factor1 * dcomplex_t{0.0, -1.0};  

        dcomplex_t factor2{0.0, factor1}; 
        for (int m = 1; m <= mmax; m += 2) {
          dcomplex_t temp1 = legendre_e[midx(n, m)] * factor2; 
          dcomplex_t temp2 = legendre_e[midx(n, m + 1)] * factor2; 
          // L_n^m where m is odd 
          W1[lidx(n, m)] += z[m] * temp1; 
          // L_n^(-m) where m is odd
          W1[lidx(n, -m)] += conj(z[m]) * temp1;
          // L_n^m where m is even
          W1[lidx(n, m + 1)] += z[m + 1] * temp2;
          // L_n^(-m) where m is even 
          W1[lidx(n, -m - 1)] += conj(z[m + 1]) * temp2;
          factor2 *= -1;
        }
      }
      delete [] z; 
    }

    // Convert propagating wave 
    for (int k = 0; k < s_p; ++k) {
      // Note: Prop(k, j) needs to be scaled by w_p[k] / m_p[k]
      int curr = 2 * sm_p[k] - 1; 
      int mk2 = m_p[k] / 2; 
          
      // Computes sum_{j = 1}^{m_p[k]} Prop(k, j) e^{-i * m * alpha_j} 
      dcomplex_t *z = new dcomplex_t[f_p[k] + 1](); 

      // m = 0
      for (int j = 1; j <= m_p[k]; ++j) {
        z[0] += Prop[curr + j]; 
      }

      // If j2 = j1 + m_p[k] / 2, then 
      // e^{-i * m * alpha_j2} = e^{-i * m * alpha_j1} * (-1)^m 
      
      // m = 1, ..., f_p[k], where m is odd 
      for (int m = 1; m <= f_p[k]; m += 2) {
        for (int j = 1; j <= mk2; ++j) {
          int aidx = smf_p[k] + (j - 1) * f_p[k] + m - 1; 
          z[m] += 
            conj(ealphaj_p[aidx]) * (Prop[curr + j] - Prop[curr + j + mk2]);
        }
      }

      // m = 2, ..., f_p[k], where m is even 
      for (int m = 2; m <= f_p[k]; m += 2) {
        for (int j = 1; j <= mk2; ++j) {
          int aidx = smf_p[k] + (j - 1) * f_p[k] + m - 1; 
          z[m] += 
            conj(ealphaj_p[aidx]) * (Prop[curr + j] + Prop[curr + j + mk2]);
        }
      }

      legendre_Plm_prop_scaled(p, cos(x_p[k]), scale, legendre_p); 
      
      double factor1 = w_p[k] / m_p[k]; 

      for (int n = 0; n <= p; ++n) {
        int mmax = (n <= f_p[k] ? n : f_p[k]); 

        // Process L_n^0
        W1[lidx(n, 0)] += z[0] * legendre_p[midx(n, 0)] * factor1; 

        for (int m = 1; m <= mmax; ++m) {
          // L_n^m
          W1[lidx(n, m)] += z[m] * legendre_p[midx(n, m)] * factor1; 
          // L_n^(-m)
          W1[lidx(n, -m)] += conj(z[m]) * legendre_p[midx(n, m)] * factor1;
        }
      }
      
      delete [] z;
    }

    // Scale the local expansion by 
    // (2 * n + 1) * (n - |m|)! / (n + |m|)! (-1)^n 
    int offset = 0; 
    for (int n = 0; n <= p; ++n) {
      double power_m1 = 1; 
      for (int m = -n; m <= n; ++m) {
        W1[offset++] *= sqf[midx(n, fabs(m))] * power_m1; 
      }
      power_m1 *= -1; 
    }

    if (!sgn) {
      // If the exponential expansion is not along the positive axis 
      // direction with respect to the source, flip the sign of the converted
      // L_n^m where n is odd 
      int offset = 1; 
      for (int n = 1; n <= p; n += 2) {
        for (int m = 0; m <= n; ++m) {
          W1[offset++] *= -1;
        }
        offset += (2 * n + 3);
      }
    }

    if (dir == 'z') {
      contrib = W1; 
    } else if (dir == 'y') {
      const double *d = builtin_helmholtz_table_->dmat_plus(0.0); 
      rotate_sph_y(W1, d, W2, true); 
      rotate_sph_z(W2, M_PI / 2, W1, true); 
      contrib = W1; 
    } else if (dir == 'x') {
      const double *d = builtin_helmholtz_table_->dmat_minus(0.0); 
      rotate_sph_y(W1, d, W2, true); 
      contrib = W2; 
    } 

    // Merge converted local expansion with the stored one
    offset = 0; 
    for (int n = 0; n <= p; ++n) {
      for (int m = -n; m <= n; ++m) {
        L[offset] += contrib[offset]; 
        offset++;
      }
    }

    delete [] W1;
    delete [] W2; 
    delete [] legendre_e; 
    delete [] legendre_p; 
  }  
}; 

} // namespace dashmm

#endif // __DASHMM_HELMHOLTZ_EXPANSION_H__
