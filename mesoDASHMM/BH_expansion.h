#ifndef __BH_EXPANSION_CPP__
#define __BH_EXPANSION_CPP__


#include <cmath>
#include <complex>
#include <memory>
#include <vector>

#include "expansion.h"
#include "index.h"
#include "particle.h"
#include "point.h"


namespace dashmm {


class BH_Expansion : public Expansion {
 public:
  BH_Expansion(Point center)
    : mtot_{0.0}, xcom_{0.0, 0.0, 0.0},
      Q_{0.0, 0.0, 0.0, 0.0, 0.0, 0.0} { }
  ~BH_Expansion() { }

  size_t size() const override {return 10;}
  Point center() const override {
    return Point{xcom_[0], xcom_[1], xcom_[2]};
  }
  std::complex<double> term(size_t i) const override;


  std::unique_ptr<Expansion> S_to_M(Point center,
      std::vector<Source>::iterator first,
      std::vector<Source>::iterator last) const override;
  std::unique_ptr<Expansion> S_to_L(Point center,
      std::vector<Source>::iterator first,
      std::vector<Source>::iterator last) const override {
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

  void M_to_T(std::vector<Target>::iterator first,
              std::vector<Target>::iterator last) const override;
  void L_to_T(std::vector<Target>::iterator first,
              std::vector<Target>::iterator last) const override {  }
  void S_to_T(std::vector<Source>::iterator s_first,
              std::vector<Source>::iterator s_last,
              std::vector<Target>::iterator t_first,
              std::vector<Target>::iterator t_last) const override;

  void add_expansion(const Expansion *temp1) override;
  void from_sum(const std::vector<const Expansion *> &exps) override;

  std::unique_ptr<Expansion> get_new_expansion(Point center) const override {
    return std::unique_ptr<Expansion>{new BH_Expansion{center}};
  }

 private:
  void calc_mtot(std::vector<Source>::iterator first,
                 std::vector<Source>::iterator last);
  void calc_xcom(std::vector<Source>::iterator first,
                 std::vector<Source>::iterator last);
  void calc_Q(std::vector<Source>::iterator first,
              std::vector<Source>::iterator last);

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
  //The ordering here is Qxx Qxy Qxz Qyy Qyz Qzz
  double Q_[6];
};


} //namespace dashmm

#endif // __BH_EXPANSION_CPP__
