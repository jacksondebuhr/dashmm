#ifndef __PARTICLE_H__
#define __PARTICLE_H__


#include "hpx/hpx.h"


//TODO likely this is a subclass of an abstract base
// Also, I am not sure how relevant this concept is. It seems useful now, but
// is it just an awkward fix for this problem, or is there something more to 
// this.
//It is a bit like a description of a possible continuation. Basically, this
// object provides an to which something can be continued.
//Another way to think about it, is as a handle that another object can use
// to call an action on another. This basically is a way to hide those HPX
// details from the user.
class SetApproxContinuation {
 public:
  SetApproxContinuation(hpx_addr_t targ, hpx_action_t act)
      : target_{targ}, action_{act} {}

  hpx_addr_t target() const {return target_;}
  hpx_action_t action() const {return action_;}
  
  //a method to explicitly continue the desired result?
  
 private:
  hpx_addr_t target_;
  hpx_action_t action_;
};


class Particle {
 public:
  Particle() : x_{0.0}, m_{0.0}, approx_{0.0}, data_{HPX_NULL} {}
  Particle(double x, double m) : x_{x}, m_{m}, approx_{0.0}, data_{HPX_NULL} {}
  
  double x() const {return x_;}
  double m() const {return m_;}
  double approx() const {return approx_;}
  
  hpx_addr_t data() const {return data_;}
  void set_data(hpx_addr_t data) {data_ = data;}
  
  void set_approx(double phi) {approx_ = phi;}
  SetApproxContinuation set_approx_cont() const;
 
 private:
  double x_;
  double m_;
  double approx_;
  hpx_addr_t data_;
};


bool operator<(const Particle &a, const Particle &b);


#endif
