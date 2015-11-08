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


typedef Expansion *(*expansion_creation_function_t)(Point);
typedef Expansion *(*expansion_interpret_function_t)(void *data, size_t bytes);


class Expansion {
 public:
  virtual ~Expansion() { }

  virtual void *release() = 0;
  virtual size_t bytes() const = 0;
  virtual bool valid() const = 0;
  virtual int type() const = 0;
  virtual bool provides_L() const = 0;
  virtual bool provides_exp() const = 0;
  virtual size_t size() const = 0;
  virtual Point center() const = 0;
  virtual std::complex<double> term(size_t i) const = 0;

  virtual std::unique_ptr<Expansion> S_to_M(Point center, Source *first,
                                            Source *last) const = 0;
  virtual std::unique_ptr<Expansion> S_to_L(Point center, Source *first,
                                            Source *last) const = 0;

  virtual std::unique_ptr<Expansion> M_to_M(int from_child,
                                            double s_size) const = 0;
  virtual std::unique_ptr<Expansion> M_to_L(Index s_index, double s_size,
                                             Index t_index) const = 0;
  virtual std::unique_ptr<Expansion> L_to_L(int to_child,
                                            double t_size) const = 0;

  virtual void M_to_T(Target *first, Target *last) const = 0;
  virtual void L_to_T(Target *first, Target *last) const = 0;
  virtual void S_to_T(Source *s_first, Source *s_last,
                      Target *t_first, Target *t_last) const = 0;

  virtual void add_expansion(const Expansion *temp1) = 0;

  virtual std::unique_ptr<Expansion> get_new_expansion(Point center) const = 0;
};


//return true on success
bool register_expansion(int type, hpx_action_t creator,
                                  hpx_action_t interpreter);


//NOTE: not intended for end-user use

//This first should be able to support an "empty" version of the object
//Or more accurately, and invalid expansion. It will have no data, but will
// be able to perform the correct methods
std::unique_ptr<Expansion> interpret_expansion(int type, void *data,
                                               size_t size);
std::unique_ptr<Expansion> create_expansion(int type, Point center);


} // namespace dashmm


#endif // __DASHMM_EXPANSION_H__
