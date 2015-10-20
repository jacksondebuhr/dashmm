#ifndef __DASHMM_EXPANSION_REF_H__
#define __DASHMM_EXPANSION_REF_H__


#include <memory>
#include <vector>

#include <hpx/hpx.h>

#include "include/expansion.h"
#include "include/index.h"
#include "include/particle.h"
#include "include/point.h"


namespace dashmm {


class ExpansionRef {
 public:
  ExpansionRef(hpx_addr_t addr) : data_{ref} { }

  bool valid() const {return data_ != HPX_NULL;}

  int type() const;

  bool provides_L() const;
  size_t size() const;
  Point center() const;

  std::complex<double> term(size_t i) const;

  std::unique_ptr<Expansion> S_to_M(Point center,
                                  std::vector<Source>::iterator first,
                                  std::vector<Source>::iterator last) const;
  std::unique_ptr<Expansion> S_to_L(Point center,
                                  std::vector<Source>::iterator first,
                                  std::vector<Source>::iterator last) const;

  std::unique_ptr<Expansion> M_to_M(int from_child, double s_size) const;
  std::unique_ptr<Expansion> M_to_L(Index s_index, double s_size,
                                    Index t_index) const;
  std::unique_ptr<Expansion> L_to_L(int to_child, double t_size) const;

  void M_to_T(std::vector<Target>::iterator first,
              std::vector<Target>::iterator last) const;
  void L_to_T(std::vector<Target>::iterator first,
              std::vector<Target>::iterator last) const;
  void S_to_T(std::vector<Source>::iterator s_first,
              std::vector<Source>::iterator s_last,
              std::vector<Target>::iterator t_first,
              std::vector<Target>::iterator t_last) const;

  void add_expansion(const Expansion *temp1);
  void from_sum(const std::vector<const Expansion *> &exps);

  std::unique_ptr<Expansion> get_new_expansion(Point center) const;

 private:
  ExpansionSerial *pin();
  void unpin();

  hpx_addr_t data_;
};


ExpansionRef globalize_expansion(Expansion *exp, hpx_addr_t where);


} // namespace dashmm


#endif // __DASHMM_EXPANSION_REF_H__
