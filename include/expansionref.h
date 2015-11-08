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
  ExpansionRef(int type, hpx_addr_t addr) : type_{type}, data_{addr} { }

  hpx_addr_t data() const {return data_;}
  bool valid() const {return data_ != HPX_NULL;}
  int type() const {return type_;}

  //NOTE: We have removed some of the normal interface for expansions...
  // I think they will not be needed, and so we have removed them.
  //And for that matter, what do we do with the rest of these? We want to
  // basically always schedule these operations. So do we remove them as well?

  //NOTE: These source and target things should be changed to SourceRef and
  // TargetRef?

  //NOTE: These do not need to wait on the expansion
  std::unique_ptr<Expansion> S_to_M(Point center,
                                    Source *first, Source *last) const;
  std::unique_ptr<Expansion> S_to_L(Point center,
                                    Source *first, Source *last) const;

  //NOTE: So do we take the set that needs to wait and make these the
  // only work through schedule
  //NOTE: These *do* have to wait for the expansion
  std::unique_ptr<Expansion> M_to_M(int from_child, double s_size) const;
  std::unique_ptr<Expansion> M_to_L(Index s_index, double s_size,
                                    Index t_index) const;
  std::unique_ptr<Expansion> L_to_L(int to_child, double t_size) const;

  //NOTE: These *do* have to wait
  void M_to_T(Target *first, Target *last) const;
  void L_to_T(Target *first, Target *last) const;

  //NOTE: This does not
  void S_to_T(Source *s_first, Source *s_last,
              Target *t_first, Target *t_last) const;

  //NOTE: This needs to wait on the input expansion
  void add_expansion(const Expansion *temp1);
  //end TODO comment

  std::unique_ptr<Expansion> get_new_expansion(Point center) const;

  //schedule is added to this

 private:
  int type_;
  hpx_addr_t data_;     //this is the LCO
};


ExpansionRef globalize_expansion(std::unique_ptr<Expansion> exp);


} // namespace dashmm


#endif // __DASHMM_EXPANSION_REF_H__
