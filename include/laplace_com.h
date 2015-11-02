#ifndef __DASHMM_LAPLACE_COM_EXPANSION_H__
#define __DASHMM_LAPLACE_COM_EXPANSION_H__


#include <cmath>
#include <complex>
#include <memory>
#include <vector>

#include "include/builtin_ids.h"
#include "include/expansion.h"
#include "include/index.h"
#include "include/particle.h"
#include "include/point.h"


namespace dashmm {


class LaplaceCOM : public Expansion {
 public:
  explicit LaplaceCOM(Point center)
      : mtot_{0.0}, xcom_{0.0, 0.0, 0.0}, Q_{0.0, 0.0, 0.0, 0.0, 0.0, 0.0} { }

  int type() const override {return kExpansionLaplaceCOM;}
  ExpansionSerialPtr serialize() const override;

  bool provides_L() const override;
  size_t size() const override {return 10;}
  Point center() const override {
    return Point{xcom_[0], xcom_[1], xcom_[2]};
  }
  std::complex<double> term(size_t i) const override;


  std::unique_ptr<Expansion> S_to_M(Point center,
                                    Source *first,
                                    Source *last) const override;
  std::unique_ptr<Expansion> S_to_L(Point center,
                                    Source *first,
                                    Source *last) const override {
    return std::unique_ptr<Expansion>{nullptr};
  }

  std::unique_ptr<Expansion> M_to_M(int from_child,
                                    double s_size) const override;
  std::unique_ptr<Expansion> M_to_L(Index s_index, double size,
                                    Index t_index) const override {
    return std::unique_ptr<Expansion>{nullptr};
  }
  std::unique_ptr<Expansion> L_to_L(int to_child,
                                    double t_size) const override {
    return std::unique_ptr<Expansion>{nullptr};
  }

  void M_to_T(Target *first, Target *last) const override;
  void L_to_T(Target *first, Target *last) const override {  }
  void S_to_T(Source *s_first, Source *s_last,
              Target *t_first, Target *t_last) const override;

  void add_expansion(const Expansion *temp1) override;
  void from_sum(const std::vector<const Expansion *> &exps) override;

  std::unique_ptr<Expansion> get_new_expansion(Point center) const override {
    return std::unique_ptr<Expansion>{new LaplaceCOM{center}};
  }

 private:
  void calc_mtot(Source *first, Source *last);
  void calc_xcom(Source *first, Source *last);
  void calc_Q(Source *first, Source *last);

  void set_mtot(double m) {mtot_ = m;}
  void set_xcom(const double c[3]) {
    xcom_[0] = c[0]; xcom_[1] = c[1]; xcom_[2] = c[2];
  }
  void set_Q(const double q[6]) {
    Q_[0] = q[0]; Q_[1] = q[1]; Q_[2] = q[2];
    Q_[3] = q[3]; Q_[4] = q[4]; Q_[5] = q[5];
  }

  double mtot_;
  double xcom_[3];
  double Q_[6];
};


} //namespace dashmm


#endif // __DASHMM_LAPLACE_COM_EXPANSION_H__
