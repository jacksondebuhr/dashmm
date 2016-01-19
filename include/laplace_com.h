#ifndef __DASHMM_LAPLACE_COM_EXPANSION_H__
#define __DASHMM_LAPLACE_COM_EXPANSION_H__


/// \file include/laplace_com.h
/// \brief Declaration of LaplaceCOM


#include <cmath>
#include <complex>
#include <memory>
#include <vector>

#include "include/expansion.h"
#include "include/ids.h"
#include "include/index.h"
#include "include/particle.h"
#include "include/point.h"


namespace dashmm {


struct LaplaceCOMData {
  int reserved;
  int type;
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
/// Largely speaking, this class implements the Expansion interface, and so
/// many of the methods are documented in Expansion.
class LaplaceCOM : public Expansion {
 public:
  explicit LaplaceCOM(Point center);
  LaplaceCOM(LaplaceCOMData *ptr, size_t bytes)
      : data_{ptr}, bytes_{sizeof(LaplaceCOMData)} { }
  ~LaplaceCOM();

  void *release() override {
    LaplaceCOMData *retval = data_;
    data_ = nullptr;
    return retval;
  }
  size_t bytes() const override {return bytes_;}
  bool valid() const override {return data_ != nullptr;}
  int type() const override {return kExpansionLaplaceCOM;}
  int accuracy() const override {return -1;} 
  bool provides_L() const override {return false;}
  bool provides_exp() const override {return false;}
  size_t size() const override {return 10;}
  Point center() const override {
    assert(valid());
    return Point{data_->xcom[0], data_->xcom[1], data_->xcom[2]};
  }
  dcomplex_t term(size_t i) const override;


  std::unique_ptr<Expansion> S_to_M(Point center, Source *first,
                                    Source *last, double scale) const override;
  std::unique_ptr<Expansion> S_to_L(Point center, Source *first,
                                    Source *last, double scale) const override 
  {
    return std::unique_ptr<Expansion>{new LaplaceCOM{nullptr, 0}};
  }

  std::unique_ptr<Expansion> M_to_M(int from_child,
                                    double s_size) const override;
  std::unique_ptr<Expansion> M_to_L(Index s_index, double size,
                                    Index t_index) const override {
    return std::unique_ptr<Expansion>{new LaplaceCOM{nullptr, 0}};
  }
  std::unique_ptr<Expansion> L_to_L(int to_child,
                                    double t_size) const override {
    return std::unique_ptr<Expansion>{new LaplaceCOM{nullptr, 0}};
  }

  void M_to_T(Target *first, Target *last, double scale) const override;
  void L_to_T(Target *first, Target *last, double scale) const override {  }
  void S_to_T(Source *s_first, Source *s_last,
              Target *t_first, Target *t_last) const override;

  void add_expansion(const Expansion *temp1) override;

  std::unique_ptr<Expansion> get_new_expansion(Point center) const override {
    return std::unique_ptr<Expansion>{new LaplaceCOM{center}};
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
    data_->xcom[0] = c[0]; data_->xcom[1] = c[1]; data_->xcom[2] = c[2];
  }

  /// Set the quadrupole moments
  ///
  /// This resets the quadrupole moments for the expansion.
  ///
  /// \param q - the new quadrupole moments
  void set_Q(const double q[6]) {
    data_->Q[0] = q[0]; data_->Q[1] = q[1]; data_->Q[2] = q[2];
    data_->Q[3] = q[3]; data_->Q[4] = q[4]; data_->Q[5] = q[5];
  }

 private:
  /// This will compute the total mass of a set of sources
  ///
  /// \param first - the first source
  /// \param last - (one past the) last source
  void calc_mtot(Source *first, Source *last);

  /// This will compute the center of mass of a set of sources
  ///
  /// \param first - the first source
  /// \param last - (one past the) last source
  void calc_xcom(Source *first, Source *last);

  /// This will compute the quadrupole moments of a set of sources
  ///
  /// \param first - the first source
  /// \param last - (one past the) last source
  void calc_Q(Source *first, Source *last);

  LaplaceCOMData *data_;
  size_t bytes_;
};


} //namespace dashmm


#endif // __DASHMM_LAPLACE_COM_EXPANSION_H__
