#ifndef __DASHMM_PARTICLE_H__
#define __DASHMM_PARTICLE_H__


#include <complex>

#include "include/point.h"


namespace dashmm {


struct Source {
  double charge;
  Point position;
};


class SourceRef {
 public:
  SourceRef() : data_{HPX_NULL}, n_{0} { }
  SourceRef(Source *sources, int n);
  SourceRef(hpx_addr_t data, int n) : data_{data}, n_{n} { }

  void destroy();
  int n() const {return n_;}
  hpx_addr_t data() const {return data_;}

 private:
  hpx_addr_t data_;
  int n_;
};


struct Target {
  Point position;
  std::complex<double> phi;
  size_t index;              //the index in the original data from the user
};


//NOTE: This will only work in SMP
class TargetRef {
 public:
  TargetRef() : data_{HPX_NULL}, n_{0} { }
  TargetRef(hpx_addr_t data, int n) : data_{data}, n_{n} { }
  TargetRef(Target *targets, int n);

  void destroy();

  int n() const {return n_;}
  hpx_addr_t data() const {return data_;}

  //TODO: some methods to schedule and so forth the various operations
  // these will likely look like the same names and stuff.
  // Perhaps these will be called by the expansionref object. This way
  // the user will not ever see this happening.
  // Or maybe the expansionref object will go ahead an do the whole thing...
  void finalize() const;
  void schedule(int num) const;

 private:
  hpx_addr_t data_;
  int n_;
};


} //namespace dashmm


#endif // __DASHMM_PARTICLE_H__
