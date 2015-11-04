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


//NOTE: Some of these will not be in the documentation.
class ExpansionRef {
 public:
  ExpansionRef(hpx_addr_t addr) : data_{ref}, exp_{nullptr} { }

  bool valid() const {return data_ != HPX_NULL;}

  int type() const;
  hpx_addr_t data() const {return data_;}

  bool provides_L() const;
  bool provides_exp() const;
  size_t size() const;
  Point center() const;

  std::complex<double> term(size_t i) const;

  std::unique_ptr<Expansion> S_to_M(Point center,
                                    Source *first, Source *last) const;
  std::unique_ptr<Expansion> S_to_L(Point center,
                                    Source *first, Source *last) const;

  std::unique_ptr<Expansion> M_to_M(int from_child, double s_size) const;
  std::unique_ptr<Expansion> M_to_L(Index s_index, double s_size,
                                    Index t_index) const;
  std::unique_ptr<Expansion> L_to_L(int to_child, double t_size) const;

  void M_to_T(Target *first, Target *last) const;
  void L_to_T(Target *first, Target *last) const;
  void S_to_T(Source *s_first, Source *s_last,
              Target *t_first, Target *t_last) const;

  void add_expansion(const Expansion *temp1);
  void from_sum(const std::vector<const Expansion *> &exps);

  std::unique_ptr<Expansion> get_new_expansion(Point center) const;

 private:
  ExpansionSerial *pin();
  void unpin();
  void setup_local_expansion();
  void save_to_global();

  hpx_addr_t data_;
  std::unique_ptr<Expansion> exp_;
};


ExpansionRef globalize_expansion(Expansion *exp, hpx_addr_t where);


} // namespace dashmm


#endif // __DASHMM_EXPANSION_REF_H__
