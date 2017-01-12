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


#ifndef __DASHMM_LAPLACE_COM_EXPANSION_H__
#define __DASHMM_LAPLACE_COM_EXPANSION_H__


/// \file
/// \brief Declaration of LaplaceCOM


#include <cassert>
#include <cmath>
#include <complex>
#include <memory>
#include <vector>

#include "dashmm/index.h"
#include "dashmm/point.h"
#include "dashmm/types.h"
#include "dashmm/viewset.h"


namespace dashmm {


struct LaplaceCOMData {
  double mtot;
  double xcom[3];
  double Q[6];
};


/// Laplace kernel Center of Mass expansion
///
/// This expansion is of the Laplace Kernel about the center of mass of the
/// represented sources. The kernel does not include any scaling for physical
/// constants, and so the user will need to multiply results of this expansion
/// by the relevant factors (including a minus sign if needed).
///
/// This expansion is most relevant for interactions that have only one sign
/// of the charge. In particular, this is especially useful for the
/// gravitational interaction. The expansion contains up to the quadrupole
/// terms. This object does not include the capability to compute local versions
/// of the expansion, so it is only compatible with methods that do not require
/// local expansions.
///
/// This class is a template with parameters for the source and target
/// types.
///
/// Source must define a double valued 'charge' member to be used with
/// LaplaceCOM. Target must define a std::complex<double> valued 'phi' member
/// to be used with LaplaceCOM.
template <typename Source, typename Target>
class LaplaceCOM {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = LaplaceCOM<Source, Target>;

  LaplaceCOM(ExpansionRole role, double scale = 1.0, Point center = Point{}) {
    bytes_ = sizeof(LaplaceCOMData);
    data_ = reinterpret_cast<LaplaceCOMData *>(new char [bytes_]);
    assert(valid(ViewSet{}));
    data_->mtot = 0.0;
    data_->xcom[0] = 0.0;
    data_->xcom[1] = 0.0;
    data_->xcom[2] = 0.0;
    data_->Q[0] = 0.0;
    data_->Q[1] = 0.0;
    data_->Q[2] = 0.0;
    data_->Q[3] = 0.0;
    data_->Q[4] = 0.0;
    data_->Q[5] = 0.0;
  }

  LaplaceCOM(const ViewSet &views) {
    assert(views.count() < 2);
    bytes_ = sizeof(LaplaceCOMData);
    if (views.count() == 1) {
      data_ = reinterpret_cast<LaplaceCOMData *>(views.view_data(0));
    } else {
      data_ = nullptr;
    }
  }

  ~LaplaceCOM() {
    if (valid(ViewSet{})) {
      delete [] data_;
      data_ = nullptr;
    }
  }

  void release() {
    data_ = nullptr;
  }

  bool valid(const ViewSet &view) const {
    assert(view.count() < 2);
    return data_ != nullptr;
  }

  int view_count() const {
    if (data_) return 1;
    return 0;
  }

  void get_views(ViewSet &view) const {
    assert(view.count() < 2);
    if (view.count() > 0) {
      view.set_bytes(0, sizeof(LaplaceCOMData));
      view.set_data(0, (char *)data_);
    }
    view.set_role(kSourcePrimary);
  }

  ViewSet get_all_views() const {
    ViewSet retval{};
    retval.add_view(0);
    get_views(retval);
    return retval;
  }

  int accuracy() const {return -1;}

  ExpansionRole role() const {return kSourcePrimary;}

  Point center() const {
    assert(valid(ViewSet{}));
    return Point{data_->xcom[0], data_->xcom[1], data_->xcom[2]};
  }

  size_t view_size(int view) const {
    assert(view == 0);
    return 10;
  }

  dcomplex_t view_term(int view, size_t i) const {
    assert(view == 0);
    if (i == 0) {
      return dcomplex_t{data_->mtot};
    } else if (i < 4) {
      return dcomplex_t{data_->xcom[i - 1]};
    } else {
      return dcomplex_t{data_->Q[i - 4]};
    }
  }

  void S_to_M(Source *first, Source *last) const {
    calc_mtot(first, last);
    calc_xcom(first, last);
    calc_Q(first, last);
  }

  std::unique_ptr<expansion_t> S_to_L(Point center, Source *first,
                                      Source *last) const {
    return std::unique_ptr<expansion_t>{nullptr};
  }

  std::unique_ptr<expansion_t> M_to_M(int from_child) const {
    assert(valid(ViewSet{}));
    expansion_t *temp = new expansion_t{kSourcePrimary}; 
    temp->set_mtot(data_->mtot);
    temp->set_xcom(data_->xcom);
    temp->set_Q(data_->Q);
    return std::unique_ptr<expansion_t>{temp};
  }

  std::unique_ptr<expansion_t> M_to_L(Index s_index, Index t_index) const {
    return std::unique_ptr<expansion_t>{nullptr};
  }

  std::unique_ptr<expansion_t> L_to_L(int to_child) const {
    return std::unique_ptr<expansion_t>{nullptr};
  }

  void M_to_T(Target *first, Target *last) const {
    assert(valid(ViewSet{}));
    for (auto i = first; i != last; ++i) {
      Point pos{i->position};

      double sum = 0;
      double diff[3] {pos.x() - data_->xcom[0], pos.y() - data_->xcom[1],
             pos.z() - data_->xcom[2]};
      double diff2mag{diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]};
      double diffmag{sqrt(diff2mag)};
      double quaddenom{-1.0 / (2.0 * diff2mag * diff2mag * diffmag)};

      sum = data_->mtot / diffmag;
      double qsum{data_->Q[0] * diff[0] * diff[0]};
      qsum += 2.0 * data_->Q[1] * diff[0] * diff[1];
      qsum += 2.0 * data_->Q[2] * diff[0] * diff[2];
      qsum += data_->Q[3] * diff[1] * diff[1];
      qsum += 2.0 * data_->Q[4] * diff[1] * diff[2];
      qsum += data_->Q[5] * diff[2] * diff[2];
      qsum *= quaddenom;
      sum += qsum;

      i->phi += dcomplex_t{sum};
    }
  }

  void L_to_T(Target *first, Target *last) const { }

  void S_to_T(Source *s_first, Source *s_last,
              Target *t_first, Target *t_last) const {
    for (auto targ = t_first; targ != t_last; ++targ) {
      Point pos = targ->position;
      double sum{0.0};
      for (auto i = s_first; i != s_last; ++i) {
        double diff[3] {pos.x() - i->position.x(),
               pos.y() - i->position.y(), pos.z() - i->position.z()};
        double mag{diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]};
        if (mag > 0) {
          sum += i->charge / sqrt(mag);
        }
      }

      targ->phi += dcomplex_t{sum};
    }
  }

  std::unique_ptr<expansion_t> M_to_I() const {
    return std::unique_ptr<expansion_t>{nullptr};
  }

  std::unique_ptr<expansion_t> I_to_I(Index s_index, Index t_index) const {
    return std::unique_ptr<expansion_t>{nullptr};
  }

  std::unique_ptr<expansion_t> I_to_L(Index t_index) const {
    return std::unique_ptr<expansion_t>{nullptr};
  }

  void add_expansion(const expansion_t *temp1) {
    double M2 = temp1->view_term(0, 0).real();
    double D2[3] = {temp1->view_term(0, 1).real(),
                    temp1->view_term(0, 2).real(),
                    temp1->view_term(0, 3).real()};
    double Q2[6] = {temp1->view_term(0, 4).real(),
                    temp1->view_term(0, 5).real(),
                    temp1->view_term(0, 6).real(),
                    temp1->view_term(0, 7).real(),
                    temp1->view_term(0, 8).real(),
                    temp1->view_term(0, 9).real()};

    double Mprime = data_->mtot + M2;

    double Dprime[3] { };
    Dprime[0] = (data_->mtot * data_->xcom[0] + M2 * D2[0]) / Mprime;
    Dprime[1] = (data_->mtot * data_->xcom[1] + M2 * D2[1]) / Mprime;
    Dprime[2] = (data_->mtot * data_->xcom[2] + M2 * D2[2]) / Mprime;

    double diff1[3] = {data_->xcom[0] - Dprime[0], data_->xcom[1] - Dprime[1],
                       data_->xcom[2] - Dprime[2]};
    double diff1mag = diff1[0] * diff1[0] + diff1[1] * diff1[1]
                      + diff1[2] * diff1[2];
    double diff2[3] = {D2[0] - Dprime[0], D2[1] - Dprime[1], D2[2] - Dprime[2]};
    double diff2mag = diff2[0] * diff2[0] + diff2[1] * diff2[1]
                      + diff2[2] * diff2[2];

    double Qprime[6] { };
    Qprime[0] = data_->mtot * (3.0 * diff1[0] * diff1[0] - diff1mag)
                + data_->Q[0];
    Qprime[1] = data_->mtot * (3.0 * diff1[0] * diff1[1]) + data_->Q[1];
    Qprime[2] = data_->mtot * (3.0 * diff1[0] * diff1[2]) + data_->Q[2];
    Qprime[3] = data_->mtot * (3.0 * diff1[1] * diff1[1] - diff1mag)
                + data_->Q[3];
    Qprime[4] = data_->mtot * (3.0 * diff1[1] * diff1[2]) + data_->Q[4];
    Qprime[5] = data_->mtot * (3.0 * diff1[2] * diff1[2] - diff1mag)
                + data_->Q[5];
    Qprime[0] += M2 * (3.0 * diff2[0] * diff2[0] - diff2mag) + Q2[0];
    Qprime[1] += M2 * (3.0 * diff2[0] * diff2[1]) + Q2[1];
    Qprime[2] += M2 * (3.0 * diff2[0] * diff2[2]) + Q2[2];
    Qprime[3] += M2 * (3.0 * diff2[1] * diff2[1] - diff2mag) + Q2[3];
    Qprime[4] += M2 * (3.0 * diff2[1] * diff2[2]) + Q2[4];
    Qprime[5] += M2 * (3.0 * diff2[2] * diff2[2] - diff2mag) + Q2[5];

    set_mtot(Mprime);
    set_xcom(Dprime);
    set_Q(Qprime);
  }

  static void update_table(int n_digits, double domain_size,
                           const std::vector<double> &kernel_params) { }

  static void delete_table() { }

  static double compute_scale(Index index) {return 1.0;}

  static int weight_estimate(Operation op,
                             Index s = Index{}, Index t = Index{}) {
    return 1;
  }

  /// Set the total mass of the expansion
  ///
  /// This sets the monopole term for the expansion.
  ///
  /// \param m - the new total mass
  void set_mtot(double m) {data_->mtot = m;}

  /// Set the expansion center
  ///
  /// This resets the center of mass of the expansion to the provided values.
  ///
  /// \param c - the new center of mass
  void set_xcom(const double c[3]) {
    data_->xcom[0] = c[0];data_->xcom[1] = c[1];data_->xcom[2] = c[2];
  }

  /// Set the quadrupole moments
  ///
  /// This resets the quadrupole moments for the expansion.
  ///
  /// \param q - the new quadrupole moments
  void set_Q(const double q[6]) {
    data_->Q[0] = q[0];data_->Q[1] = q[1];data_->Q[2] = q[2];
    data_->Q[3] = q[3];data_->Q[4] = q[4];data_->Q[5] = q[5];
  }

 private:
  /// This will compute the total mass of a set of sources
  ///
  /// \param first - the first source
  /// \param last - (one past the) last source
  void calc_mtot(Source *first, Source *last) const {
    assert(valid(ViewSet{}));
    data_->mtot = 0.0;
    for (auto i = first; i != last; ++i) {
      data_->mtot += i->charge;
    }
  }

  /// This will compute the center of mass of a set of sources
  ///
  /// \param first - the first source
  /// \param last - (one past the) last source
  void calc_xcom(Source *first, Source *last) const {
    assert(valid(ViewSet{}));
    data_->xcom[0] = 0.0;
    data_->xcom[1] = 0.0;
    data_->xcom[2] = 0.0;
    for (auto i = first; i != last; ++i) {
      data_->xcom[0] += i->charge * i->position.x();
      data_->xcom[1] += i->charge * i->position.y();
      data_->xcom[2] += i->charge * i->position.z();
    }
    if (data_->mtot != 0.0) {
      double oomtot = 1.0 / data_->mtot;
      data_->xcom[0] *= oomtot;
      data_->xcom[1] *= oomtot;
      data_->xcom[2] *= oomtot;
    }
  }

  /// This will compute the quadrupole moments of a set of sources
  ///
  /// \param first - the first source
  /// \param last - (one past the) last source
  void calc_Q(Source *first, Source *last) const {
    assert(valid(ViewSet{}));
    for (int i = 0; i < 6; ++i) {
      data_->Q[i] = 0;
    }

    for (auto i = first; i != last; ++i) {
      double diff[3] = {i->position.x() - data_->xcom[0],
                        i->position.y() - data_->xcom[1],
                        i->position.z() - data_->xcom[2]};
      double diffmag{diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]};

      data_->Q[0] += i->charge * (3.0 * diff[0] * diff[0] - diffmag);   // Qxx
      data_->Q[1] += i->charge * 3.0 * diff[0] * diff[1];               // Qxy
      data_->Q[2] += i->charge * 3.0 * diff[0] * diff[2];               // Qxz
      data_->Q[3] += i->charge * (3.0 * diff[1] * diff[1] - diffmag);   // Qyy
      data_->Q[4] += i->charge * 3.0 * diff[1] * diff[2];               // Qyz
      data_->Q[5] += i->charge * (3.0 * diff[2] * diff[2] - diffmag);   // Qzz
    }
  }

  LaplaceCOMData *data_;
  size_t bytes_;
};


} // namespace dashmm


#endif // __DASHMM_LAPLACE_COM_EXPANSION_H__
