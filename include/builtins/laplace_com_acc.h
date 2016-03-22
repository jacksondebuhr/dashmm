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


#ifndef __DASHMM_LAPLACE_COM_ACC_EXPANSION_H__
#define __DASHMM_LAPLACE_COM_ACC_EXPANSION_H__


/// \file include/laplace_com_acc.h
/// \brief Declaration of LaplaceCOMAcc


#include <cmath>
#include <complex>
#include <memory>
#include <vector>

#include "dashmm/index.h"
#include "dashmm/point.h"
#include "dashmm/types.h"


// NOTE: Built-in Methods and Expansions will be part of the dashmm namespace.
// User-defined Methods and Expansion will not be.
namespace dashmm {


struct LaplaceCOMAccData {
  int reserved;
  int n_digits; // unused
  double mtot;
  double xcom[3];
  double Q[6];
};


/// Laplace kernel Center of Mass expansion yielding acceleration
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
/// LaplaceCOMAcc. Target must define a double [3] valued 'acceleration' member
/// to be used with LaplaceCOMAcc.
template <typename Source, typename Target>
class LaplaceCOMAcc {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = LaplaceCOMAcc<Source, Target>;

  LaplaceCOMAcc(Point center, int n_digits) {
    bytes_ = sizeof(LaplaceCOMAccData);
    data_ = reinterpret_cast<LaplaceCOMAccData *>(new char [bytes_]);
    assert(valid());
    data_->n_digits = -1; // unused
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

  LaplaceCOMAcc(void *ptr, size_t bytes, int n_digits)
      : data_{static_cast<LaplaceCOMAccData *>(ptr)},
        bytes_{sizeof(LaplaceCOMAccData)} { }

  ~LaplaceCOMAcc() {
    if (valid()) {
      delete [] data_;
      data_ = nullptr;
    }
  }

  void *release() {
    LaplaceCOMAccData *retval = data_;
    data_ = nullptr;
    return retval;
  }

  size_t bytes() const {return bytes_;}

  bool valid() const {return data_ != nullptr;}

  int accuracy() const {return -1;}

  size_t size() const {return 10;}

  Point center() const {
    assert(valid());
    return Point{data_->xcom[0], data_->xcom[1], data_->xcom[2]};
  }

  dcomplex_t term(size_t i) const {
    if (i == 0) {
      return dcomplex_t{data_->mtot};
    } else if (i < 4) {
      return dcomplex_t{data_->xcom[i - 1]};
    } else {
      return dcomplex_t{data_->Q[i - 4]};
    }
  }

  void S_to_M(Point center, Source *first, Source *last,
              double scale) const {
    calc_mtot(first, last);
    calc_xcom(first, last);
    calc_Q(first, last);
  }

  std::unique_ptr<expansion_t> S_to_L(Point center,
                                      Source *first, Source *last,
                                      double scale) const {
    return std::unique_ptr<expansion_t>{nullptr};
  }

  std::unique_ptr<expansion_t> M_to_M(int from_child,
                                      double s_size) const {
    assert(valid());
    expansion_t *temp = new expansion_t(Point{0.0, 0.0, 0.0}, 0);
    temp->set_mtot(data_->mtot);
    temp->set_xcom(data_->xcom);
    temp->set_Q(data_->Q);
    return std::unique_ptr<expansion_t>{temp};
  }

  std::unique_ptr<expansion_t> M_to_L(Index s_index, double s_size,
                                      Index t_index) const {
    return std::unique_ptr<expansion_t>{nullptr};
  }

  std::unique_ptr<expansion_t> L_to_L(int to_child,
                                      double t_size) const {
    return std::unique_ptr<expansion_t>{nullptr};
  }

  void M_to_T(Target *first, Target *last, double scale) const {
    assert(valid());
    // TODO
    for (auto i = first; i != last; ++i) {
      Point pos{i->position};

      double diff[3] = {pos.x() - data_->xcom[0], pos.y() - data_->xcom[1],
                        pos.z() - data_->xcom[2]};
      double diff2mag{diff[0] * diff[0] + diff[1] * diff[1]
                      + diff[2] * diff[2]};
      double diffmag{sqrt(diff2mag)};
      double nhat[3] = {diff[0] / diffmag, diff[1] / diffmag,
                        diff[2] / diffmag};
      double diff4mag{diff2mag * diff2mag};

      // Monopole
      double mono = data_->mtot / diff2mag;
      i->acceleration[0] += mono * nhat[0];
      i->acceleration[1] += mono * nhat[1];
      i->acceleration[2] += mono * nhat[2];

      // One quadrupole term
      double qsum{data_->Q[0] * nhat[0] * nhat[0]};
      qsum += 2.0 * data_->Q[1] * nhat[0] * nhat[1];
      qsum += 2.0 * data_->Q[2] * nhat[0] * nhat[2];
      qsum += data_->Q[3] * nhat[1] * nhat[1];
      qsum += 2.0 * data_->Q[4] * nhat[1] * nhat[2];
      qsum += data_->Q[5] * nhat[2] * nhat[2];
      qsum *= 5.0 / (2.0 * diff4mag);
      i->acceleration[0] += qsum * nhat[0];
      i->acceleration[1] += qsum * nhat[1];
      i->acceleration[2] += qsum * nhat[2];

      // The other quadrupole term
      double qalt{0.0};
      qalt = data_->Q[0] * nhat[0] + data_->Q[1] * nhat[1]
             + data_->Q[2] * nhat[2];
      i->acceleration[0] -= qalt / diff4mag;

      qalt = data_->Q[1] * nhat[0] + data_->Q[3] * nhat[1]
             + data_->Q[4] * nhat[2];
      i->acceleration[1] -= qalt / diff4mag;

      qalt = data_->Q[2] * nhat[0] + data_->Q[4] * nhat[1]
             + data_->Q[5] * nhat[2];
      i->acceleration[2] -= qalt / diff4mag;
    }
  }

  void L_to_T(Target *first, Target *last, double scale) const { }

  void S_to_T(Source *s_first, Source *s_last,
              Target *t_first, Target *t_last) const {
    for (auto targ = t_first; targ != t_last; ++targ) {
      Point pos = targ->position;
      double sum[3] = {0.0, 0.0, 0.0};
      for (auto i = s_first; i != s_last; ++i) {
        double diff[3] {pos.x() - i->position.x(),
               pos.y() - i->position.y(), pos.z() - i->position.z()};
        double mag{diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]};
        double sqrtmag{sqrt(mag)};
        if (mag > 0) {
          sum[0] += i->charge * diff[0] / (mag * sqrtmag);
          sum[1] += i->charge * diff[1] / (mag * sqrtmag);
          sum[2] += i->charge * diff[2] / (mag * sqrtmag);
        }
      }

      targ->acceleration[0] += sum[0];
      targ->acceleration[1] += sum[1];
      targ->acceleration[2] += sum[2];
    }
  }

  void add_expansion(const expansion_t *temp1) {
    double M2 = temp1->term(0).real();
    double D2[3] = {temp1->term(1).real(), temp1->term(2).real(),
                    temp1->term(3).real()};
    double Q2[6] = {temp1->term(4).real(), temp1->term(5).real(),
                    temp1->term(6).real(), temp1->term(7).real(),
                    temp1->term(8).real(), temp1->term(9).real()};

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
    data_->xcom[0] = c[0];
    data_->xcom[1] = c[1];
    data_->xcom[2] = c[2];
  }

  /// Set the quadrupole moments
  ///
  /// This resets the quadrupole moments for the expansion.
  ///
  /// \param q - the new quadrupole moments
  void set_Q(const double q[6]) {
    data_->Q[0] = q[0];
    data_->Q[1] = q[1];
    data_->Q[2] = q[2];
    data_->Q[3] = q[3];
    data_->Q[4] = q[4];
    data_->Q[5] = q[5];
  }

 private:
  /// This will compute the total mass of a set of sources
  ///
  /// \param first - the first source
  /// \param last - (one past the) last source
  void calc_mtot(Source *first, Source *last) const {
    assert(valid());
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
    assert(valid());
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
    assert(valid());
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

  LaplaceCOMAccData *data_;
  size_t bytes_;
};


} // namespace dashmm


#endif // __DASHMM_LAPLACE_COM_ACC_EXPANSION_H__
