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
#include "builtins/merge_shift.h"
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

void yuk_rotate_sph_z(const dcomplex_t *M, double alpha, dcomplex_t *MR); 
void yuk_rotate_sph_y(const dcomplex_t *M, const double *d, dcomplex_t *MR); 

void yuk_s_to_m(Point dist, double q, double scale, dcomplex_t *M);
void yuk_s_to_l(Point dist, double q, double scale, dcomplex_t *L); 
void yuk_m_to_m(int from_child, const dcomplex_t *M, double scale, 
                dcomplex_t *W); 
void yuk_l_to_l(int to_child, const dcomplex_t *L, double scale, 
                dcomplex_t *W); 
void yuk_m_to_i(const dcomplex_t *M, ViewSet &views, double scale, int id); 
void yuk_i_to_i(Index s_index, Index t_index, const ViewSet &s_views, 
                int sid, int tid, double scale, ViewSet &t_views); 
void yuk_i_to_l(const ViewSet &views, int id, Index t_index, double scale, 
                dcomplex_t *L); 
void yuk_e_to_e(dcomplex_t *M, const dcomplex_t *W, int x, int y, int z, 
                double scale); 
void yuk_e_to_l(const dcomplex_t *E, char dir, bool sgn, double scale, 
                dcomplex_t *L); 

std::vector<double> yuk_m_to_t(Point dist, double scale, 
                               const dcomplex_t *M, bool g = false);
std::vector<double> yuk_l_to_t(Point dist, double scale, 
                               const dcomplex_t *L, bool g = false);

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

  Yukawa(ExpansionRole role, double scale = 1.0, Point center = Point{})
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

  std::unique_ptr<expansion_t> S_to_M(Source *first, Source *last) const {
    double scale = views_.scale();
    Point center = views_.center();
    expansion_t *retval{new expansion_t{kSourcePrimary, scale, center}};
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, center); 
      yuk_s_to_m(dist, i->charge, scale, M); 
    }
    return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> S_to_L(Source *first, Source *last) const {
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
    double scale = views_.scale();
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    dcomplex_t *W = reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    yuk_m_to_m(from_child, M, scale, W); 
    return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> M_to_L(Index s_index, Index t_index) const {
    return std::unique_ptr<expansion_t>{nullptr};
  }

  std::unique_ptr<expansion_t> L_to_L(int to_child) const {
    expansion_t *retval{new expansion_t{kTargetPrimary}};
    double scale = views_.scale();
    dcomplex_t *L = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    dcomplex_t *W = reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    yuk_l_to_l(to_child, L, scale, W); 
    return std::unique_ptr<expansion_t>{retval};
  }

  void M_to_T(Target *first, Target *last) const {
    double scale = views_.scale();
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(views_.view_data(0));

    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, views_.center()); 
      auto result = yuk_m_to_t(dist, scale, M); 
      i->phi += result[0]; 
    }
  }

  void L_to_T(Target *first, Target *last) const {
    double scale = views_.scale();
    dcomplex_t *L = reinterpret_cast<dcomplex_t *>(views_.view_data(0));

    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, views_.center()); 
      auto result = yuk_l_to_t(dist, scale, L); 
      i->phi += result[0];
    }
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
    expansion_t *retval{new expansion_t{kSourceIntermediate, scale}};
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    yuk_m_to_i(M, retval->views_, scale, 0); 
    return std::unique_ptr<expansion_t>(retval);
  }

  std::unique_ptr<expansion_t> I_to_I(Index s_index, Index t_index) const {
    ViewSet views{kTargetIntermediate};
    double scale = views_.scale(); 
    yuk_i_to_i(s_index, t_index, views_, 0, 0, scale, views); 
    expansion_t *retval = new expansion_t{views};
    return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> I_to_L(Index t_index) const {
    // t_index is the index of the child
    expansion_t *retval{new expansion_t{kTargetPrimary}};
    dcomplex_t *L = reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    double scale = views_.scale() / 2.0; 
    yuk_i_to_l(views_, 0, t_index, scale, L); 
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
};


} // namespace dashmm

#endif // __DASHMM_YUKAWA_EXPANSION_H__




