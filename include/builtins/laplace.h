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


#ifndef __DASHMM_LAPLACE_EXPANSION_H__
#define __DASHMM_LAPLACE_EXPANSION_H__


/// \file
/// \brief Declaration of Laplace


#include <cassert>
#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <vector>

#include "dashmm/index.h"
#include "builtins/laplace_table.h"
#include "builtins/merge_shift.h"
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

void lap_rotate_sph_z(const dcomplex_t *M, double alpha, dcomplex_t *MR);
void lap_rotate_sph_y(const dcomplex_t *M, const double *d, dcomplex_t *MR);

void lap_s_to_m(Point dist, double q, double scale, dcomplex_t *M);
void lap_s_to_l(Point dist, double q, double scale, dcomplex_t *L);
void lap_m_to_m(int from_child, const dcomplex_t *M, dcomplex_t *W);
void lap_l_to_l(int to_child, const dcomplex_t *L, dcomplex_t *W);
void lap_m_to_i(const dcomplex_t *M, ViewSet &views, int id);
void lap_i_to_i(Index s_index, Index t_index, const ViewSet &s_views,
                int sid, int tid, ViewSet &t_views);
void lap_i_to_l(const ViewSet &views, int id, Index t_index, double scale,
                dcomplex_t *L);
void lap_e_to_e(dcomplex_t *M, const dcomplex_t *W, int x, int y, int z);
void lap_e_to_l(const dcomplex_t *E, char dir, bool sgn, dcomplex_t *L);

std::vector<double> lap_m_to_t(Point dist, double scale,
                               const dcomplex_t *M, bool g = false);
std::vector<double> lap_l_to_t(Point dist, double scale,
                               const dcomplex_t *L, bool g = false);


/// This class is a template with parameters for the source and target
/// types.
///
/// Source must define a double valued 'charge' member to be used with
/// Laplace. Target must define a std::complex<double> valued 'phi' member
/// to be used with Laplace.
template <typename Source, typename Target>
class Laplace {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Laplace<Source, Target>;

  Laplace(ExpansionRole role, double scale = 1.0, Point center = Point{})
    : views_{ViewSet{role, center, scale}} {
    // View size for each spherical harmonic expansion
    int p = builtin_laplace_table_->p();
    int nsh = (p + 1) * (p + 2) / 2;

    // View size for each exponential expansion
    int nexp = builtin_laplace_table_->nexp();

    if (role == kSourcePrimary || role == kTargetPrimary) {
      size_t bytes = sizeof(dcomplex_t) * nsh;
      char *data = new char[bytes]();
      views_.add_view(0, bytes, data);
    } else if (role == kSourceIntermediate) {
      size_t bytes = sizeof(dcomplex_t) * nexp;
      for (int i = 0; i < 6; ++i) {
        char *data = new char[bytes]();
        views_.add_view(i, bytes, data);
      }
    } else if (role == kTargetIntermediate) {
      size_t bytes = sizeof(dcomplex_t) * nexp;
      for (int i = 0; i < 28; ++i) {
        char *data = new char[bytes]();
        views_.add_view(i, bytes, data);
      }
    }
  }

  Laplace(const ViewSet &views) : views_{views} { }

  ~Laplace() {
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

  std::unique_ptr<expansion_t> S_to_M(const Source *first,
                                      const Source *last) const {
    double scale = views_.scale();
    Point center = views_.center();
    expansion_t *retval{new expansion_t{kSourcePrimary, scale, center}};
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, center);
      lap_s_to_m(dist, i->charge, scale, M);
    }
   return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> S_to_L(const Source *first,
                                      const Source *last) const {
    double scale = views_.scale();
    Point center = views_.center();
    expansion_t *retval{new expansion_t{kTargetPrimary}};
    dcomplex_t *L = reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, center);
      lap_s_to_l(dist, i->charge, scale, L);
    }
    return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> M_to_M(int from_child) const {
    expansion_t *retval{new expansion_t{kSourcePrimary}};
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    dcomplex_t *W = reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    lap_m_to_m(from_child, M, W);
    return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> M_to_L(Index s_index, Index t_index) const {
    expansion_t *retval{new expansion_t{kTargetPrimary}};

    int t2s_x = s_index.x() - t_index.x();
    int t2s_y = s_index.y() - t_index.y();
    int t2s_z = s_index.z() - t_index.z();
    int p = builtin_laplace_table_->p();

    // Shifting distance
    double rho = sqrt(t2s_x * t2s_x + t2s_y * t2s_y + t2s_z * t2s_z);

    // Compute powers of rho
    double *powers_rho = new double[p * 2 + 1];
    powers_rho[0] = 1.0 / rho;
    for (int i = 1; i <= p * 2; i++) {
      powers_rho[i] = powers_rho[i - 1] / rho;
    }

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
        M_to_L_zp(M, powers_rho, W1);
      } else {
        M_to_L_zm(M, powers_rho, W1);
      }
    } else {
      // azimuthal angle
      double beta = acos(t2s_x / proj);
      if (t2s_y < 0) {
        beta = 2 * M_PI - beta;
      }

      // Get precomputed Wigner d-matrix for rotation about y-axis
      const double *d1 = builtin_laplace_table_->dmat_plus(t2s_z / rho);
      const double *d2 = builtin_laplace_table_->dmat_minus(t2s_z / rho);
      lap_rotate_sph_z(M, beta, W1);
      lap_rotate_sph_y(W1, d1, W2);
      M_to_L_zp(W2, powers_rho, W1);
      lap_rotate_sph_y(W1, d2, W2);
      lap_rotate_sph_z(W2, -beta, W1);
    }

    delete [] W2;
    delete [] powers_rho;
    return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> L_to_L(int to_child) const {
    expansion_t *retval{new expansion_t{kTargetPrimary}};
    dcomplex_t *L = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    dcomplex_t *W = reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    lap_l_to_l(to_child, L, W);
    return std::unique_ptr<expansion_t>{retval};
  }

  void M_to_T(Target *first, Target *last) const {
    double scale = views_.scale();
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(views_.view_data(0));

    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, views_.center());
      auto result = lap_m_to_t(dist, scale, M);
      i->phi += result[0];
    }
  }

  void L_to_T(Target *first, Target *last) const {
    double scale = views_.scale();
    dcomplex_t *L = reinterpret_cast<dcomplex_t *>(views_.view_data(0));

    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, views_.center());
      auto result = lap_l_to_t(dist, scale, L);
      i->phi += result[0];
    }
  }

  void S_to_T(const Source *s_first,
              const Source *s_last,
              Target *t_first,
              Target *t_last) const {
    for (auto i = t_first; i != t_last; ++i) {
      dcomplex_t potential{0.0, 0.0};
      for (auto j = s_first; j != s_last; ++j) {
        Point s2t = point_sub(i->position, j->position);
        double dist = s2t.norm();
        if (dist > 0) {
          potential += j->charge / dist;
        }
      }
      i->phi += potential;
    }
  }

  std::unique_ptr<expansion_t> M_to_I() const {
    expansion_t *retval{new expansion_t{kSourceIntermediate}};
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    lap_m_to_i(M, retval->views_, 0);
    return std::unique_ptr<expansion_t>(retval);
  }

  std::unique_ptr<expansion_t> I_to_I(Index s_index, Index t_index) const {
    ViewSet views{kTargetIntermediate};
    lap_i_to_i(s_index, t_index, views_, 0, 0, views);
    expansion_t *retval = new expansion_t{views};
    return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> I_to_L(Index t_index) const {
    // t_index is the index of the child
    expansion_t *retval{new expansion_t{kTargetPrimary}};
    dcomplex_t *L = reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    double scale = views_.scale() * 2;
    lap_i_to_l(views_, 0, t_index, scale, L);
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
    update_laplace_table(n_digits, domain_size);
  }

  static void delete_table() { }

  static double compute_scale(Index index) {
    return builtin_laplace_table_->scale(index.level());
  }

  static int weight_estimate(Operation op,
                             Index s = Index{}, Index t = Index{}) {
    int weight = 0;
    if (op == Operation::MtoI) {
      weight = 6;
    } else if (op == Operation::ItoI) {
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
  void M_to_L_zp(const dcomplex_t *M, const double *rho, dcomplex_t *L) const {
    int p = builtin_laplace_table_->p();
    const double *sqbinom = builtin_laplace_table_->sqbinom();
    double scale = views_.scale();

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

  void M_to_L_zm(const dcomplex_t *M, const double *rho, dcomplex_t *L) const {
    int p = builtin_laplace_table_->p();
    const double *sqbinom = builtin_laplace_table_->sqbinom();
    double scale = views_.scale();

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


} // namespace dashmm

#endif // __DASHMM_LAPLACE_EXPANSION_H__
