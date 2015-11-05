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


//NOTE: This will only work in SMP
class SourceRef {
 public:
  SourceRef(hpx_addr_t data, int n, int n_total)
      : local_{nullptr}, data_{data}, n_{n}, n_total_{n_total} { }
  ~SourceRef() {unpin();}

  Source *first() const {
    pin();
    return local_;
  }
  Source *last() const {
    pin();
    return &local_[n_];
  }

  int n() const {return n_;}
  int n_total() const {return n_total_;}
  hpx_addr_t data() const {return data_;}

 private:
  void pin() const {
    if (!local_ && data_ != HPX_NULL) {
      assert(hpx_gas_try_pin(data_, (void **)&local_));
    }
  }
  void unpin() const {
    if (local_ && data_ != HPX_NULL) {
      hpx_gas_unpin(data_);
    }
  }

  mutable Source *local_;

  hpx_addr_t data_;
  int n_;
  int n_total_;
};


class Target {
 public:
  Target() : position_{0.0, 0.0, 0.0}, phi_{0.0}, index_{0} { }
  Target(double x, double y, double z, size_t idx) :
    position_{x, y, z}, phi_{}, direct_{}, index_{idx} { }

  Point position() const {return position_;}
  double x() const {return position_.x();}
  double y() const {return position_.y();}
  double z() const {return position_.z();}
  std::complex<double> phi() const {return phi_;}
  size_t index() const {return index_;}

  void set_position(Point p) {position_ = p;}
  void set_phi(std::complex<double> p) {phi_ = p;}
  void set_index(size_t idx) {index_ = idx;}

 private:
  Point position_;
  std::complex<double> phi_;
  size_t index_;              //the index in the original data from the user
};


//NOTE: This will only work in SMP
class TargetRef {
 public:
  TargetRef(hpx_addr_t data, int n, int n_total)
      : local_{nullptr}, data_{data}, n_{n}, n_total_{n_total} { }
  ~TargetRef() {unpin();}

  Target *first() const {
    pin();
    return local_;
  }
  Target *last() const {
    pin();
    return &local_[n_];
  }

  int n() const {return n_;}
  int n_total() const {return n_total_;}
  hpx_addr_t data() const {return data_;}

 private:
  void pin() const {
    if (!local_ && data_ != HPX_NULL) {
      assert(hpx_gas_try_pin(data_, (void **)&local_));
    }
  }
  void unpin() const {
    if (local_ && data_ != HPX_NULL) {
      hpx_gas_unpin(data_);
    }
  }

  mutable Target *local_;

  hpx_addr_t data_;
  int n_;
  int n_total_;
};


} //namespace dashmm


#endif // __DASHMM_PARTICLE_H__
