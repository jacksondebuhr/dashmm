#ifndef __DASHMM_PARTICLE_H__
#define __DASHMM_PARTICLE_H__


#include <complex>

#include "include/point.h"


namespace dashmm {


class Source {
 public:
  Source() : charge_{0.0}, position_{0.0, 0.0, 0.0} { }
  Source(double x, double y, double z, double q) :
    charge_{q}, position_{x, y, z} { }

  double charge() const {return charge_;}
  Point position() const {return position_;}
  double x() const {return position_.x();}
  double y() const {return position_.y();}
  double z() const {return position_.z();}

  void set_charge(double q) {charge_ = q;}
  void set_position(Point p) {position_ = p;}
 private:
  double charge_;
  Point position_;
};


class Target {
 public:
  Target() : position_{0.0, 0.0, 0.0}, phi_{0.0}, direct_{0.0} { }
  Target(double x, double y, double z) :
    position_{x, y, z}, phi_{}, direct_{} { }

  Point position() const {return position_;}
  double x() const {return position_.x();}
  double y() const {return position_.y();}
  double z() const {return position_.z();}
  std::complex<double> phi() const {return phi_;}
  std::complex<double> direct() const {return direct_;}

  void set_phi(std::complex<double> p) {phi_ = p;}
  void set_direct(std::complex<double> p) {direct_ = p;}

 private:
  Point position_;
  std::complex<double> phi_;
  std::complex<double> direct_;
};


} //namespace dashmm


#endif // __DASHMM_PARTICLE_H__
