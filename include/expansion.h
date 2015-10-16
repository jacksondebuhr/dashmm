#ifndef __DASHMM_EXPANSION_H__
#define __DASHMM_EXPANSION_H__


#include <cmath>
#include <complex>
#include <memory>
#include <vector>

#include <hpx/hpx.h>

#include "include/index.h"
#include "include/particle.h"


namespace dashmm {


extern constexpr int kFirstUserExpansionType;
extern constexpr int kLastUserExpansionType;


typedef Expansion *(*expansion_creation_function_t)(size_t, void *);


struct ExpansionSerial {
  int type;
  //TODO is it useful to put common expansion stuff here?
  //that way we could get just this stuff, and not all of the data if needed?
  size_t size;
  char data[];
}


using ExpansionSerialPtr =
          std::unique_ptr<ExpansionSerial, void (*)(ExpansionSerial *)>;


class Expansion {
 public:
  virtual ~Expansion() { }

  virtual int type() const = 0;
  virtual ExpansionSerialPtr serialize() const = 0;

  virtual bool provides_L() const = 0;
  virtual size_t size() const = 0;
  virtual Point center() const = 0;

  virtual std::complex<double> term(size_t i) const = 0;

  virtual std::unique_ptr<Expansion> S_to_M(Point center,
                                  std::vector<Source>::iterator first,
                                  std::vector<Source>::iterator last) const = 0;
  virtual std::unique_ptr<Expansion> S_to_L(Point center,
                                  std::vector<Source>::iterator first,
                                  std::vector<Source>::iterator last) const = 0;

  virtual std::unique_ptr<Expansion> M_to_M(int from_child,
                                            double s_size) const = 0;
  virtual std::unique_ptr<Expansion> M_to_L(Index s_index, double s_size,
                                            Index t_index) const = 0;
  virtual std::unique_ptr<Expansion> L_to_L(int to_child,
                                            double t_size) const = 0;

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


//TODO: put in something for the ExpansionRef, either the definition here
// or an inclusion of another file...


//return true on success
//TODO hpx is exposed here...
bool register_expansion(int type, hpx_action_t creator);

//this should be used by the serialize methods for expansions
ExpansionSerialPtr expansion_serialization_allocator(size_t size);

//TODO exposure of HPX here...
ExpansionRef globalize_expansion(Expansion *exp, hpx_addr_t where);


} // namespace dashmm


#endif // __DASHMM_EXPANSION_H__
