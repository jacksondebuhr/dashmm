#ifndef __EXPANSION_H__
#define __EXPANSION_H__


#include <cmath>
#include <complex>
#include <memory>
#include <vector>

#include "index.h"
#include "particle.h"


namespace dashmm {


//TODO:
// There is another opportunity for abstraction here. We have defined some of
// these methods to return a single double. But I can imagine cases where we
// produce multiple values. So there is some notion here that we can define
// to be the actual values produced. Perhaps it is sufficient to support a
// number of doubles. But perhaps there are other things that we want to
// be able to produce with DASHMM.
//

class Expansion {
 public:
  //NOTE: Construction of derived classes must look like the following.
  //Expansion(Point center);
  virtual ~Expansion() { }

  virtual size_t size() const = 0;
  virtual Point center() const = 0;

  virtual std::complex<double> term(size_t i) const = 0;

  virtual std::unique_ptr<Expansion> 
  S_to_M(Point center, std::vector<Source>::iterator first,
         std::vector<Source>::iterator last) const = 0;

  virtual std::unique_ptr<Expansion> 
  S_to_L(Point center, std::vector<Source>::iterator first,
         std::vector<Source>::iterator last) const = 0;

  virtual std::unique_ptr<Expansion> 
  M_to_M(int from_child, double s_size) const = 0;

  virtual std::unique_ptr<Expansion> 
  M_to_L(Index s_index, double s_size, Index t_index) const = 0;

  virtual std::unique_ptr<Expansion> 
  L_to_L(int to_child, double t_size) const = 0;

  virtual void M_to_T(std::vector<Target>::iterator first,
                      std::vector<Target>::iterator last) const = 0;
  virtual void L_to_T(std::vector<Target>::iterator first,
                      std::vector<Target>::iterator last) const = 0;
  virtual void S_to_T(std::vector<Source>::iterator s_first,
                      std::vector<Source>::iterator s_last,
                      std::vector<Target>::iterator t_first,
                      std::vector<Target>::iterator t_last) const = 0;

  virtual void add_expansion(const Expansion *temp1) = 0;
  virtual void from_sum(const std::vector<const Expansion *> &exps) = 0;

  virtual std::unique_ptr<Expansion> get_new_expansion(Point center) const = 0;
};


} //namespace dashmm


#endif // __EXPANSION_H__
