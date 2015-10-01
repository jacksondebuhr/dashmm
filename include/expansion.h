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


struct ExpansionDesc {
  bool provides_L;
  size_t size;
  hpx_action_t destroy_function;
  hpx_action_t S_to_M_function;
  hpx_action_t S_to_L_function;
  hpx_action_t M_to_M_function;
  hpx_action_t M_to_L_function;
  hpx_action_t L_to_L_function;
  hpx_action_t M_to_T_function;
  hpx_action_t L_to_T_function;
  hpx_action_t S_to_T_function;
  hpx_action_t add_expansion_function;
  hpx_action_t from_sum_function;
  hpx_action_t get_new_expansion_function;
  //The following is a pointer to shared data for each instance of the
  // expansion. We need a copy of this data on each locality. It will be part
  // of the lookup table about the expansion. So for expansions when we make
  // new ones, we shall always have to look up the core pointer?
  // Or maybe we can keep both the size and the pointer. If the size is nonzero,
  // then look it up.
  size_t core_size;
  void *core;
};


class Expansion {
 public:
  //NOTE: This is again generally not used by the user. The preference is that
  // user code will copy construct these objects.
  //allocate controls if the object will allocate global memory
  Expansion(int type, Point center, bool allocate = true);

  hpx_addr_t data() const {return data_;}
  bool valid() const {return data_ != HPX_NULL;}
  bool provides_L() const {return table_.provides_L;}
  size_t size() const {return table_.size;}
  int type() const {return type_;}

  void destroy();

  Point center() const;

  //TODO: Finish going over these; make sure the argument and return types are
  // solid.

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
  int type_;
  const ExpansionDesc &table_;
  hpx_addr_t data_;
};


//The function signature for coredata generating functions
typedef void *(*coregen_function_t)(size_t);

int register_expansion(ExpansionDesc desc, hpx_action_t coregen);


} // namespace dashmm


#endif // __DASHMM_EXPANSION_H__
