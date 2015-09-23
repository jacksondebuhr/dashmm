#ifndef __DASHMM_EXPANSION_H__
#define __DASHMM_EXPANSION_H__


#include <cmath>
#include <complex>
#include <memory>
#include <vector>

#include "hpx/hpx.h"

#include "include/index.h"
#include "include/particle.h"


namespace dashmm {


class Expansion {
 public:
  //NOTE: Construction of derived classes must provide the following
  //explicit Expansion(Point center);
  //explicit Expansion(hpx_addr_t data);

  virtual void destroy() = 0;

  virtual hpx_addr_t data() const = 0;
  virtual bool valid() const = 0;

  virtual bool provides_L() const = 0;

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


} // namespace dashmm


#endif // __DASHMM_EXPANSION_H__
