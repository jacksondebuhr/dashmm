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

  void destroy();

  hpx_addr_t data() const {return data_;}
  bool valid() const {return data_ != HPX_NULL;}
  int type() const {return type_;}

  //TODO this one provides the serial data for the object
  // NOTE: These will be synchronous... Only call if the expansion is ready
  void *serial() const;
  size_t bytes() const;

  //TODO: these...
  // So all of these basically have the following pattern: once the LCO is
  // ready, an action will spawn to continue the data to a relevant other
  // action.
  //

  //NOTE: These do not need to wait on the expansion
  // we will have access to the SourceRef, which knows counts. SO here we
  // just get the source data, compute the S_to_M and then set the LCO
  //
  // TODO: This needs more thinking too...
  //
  // See next. This is generally called on the prototype, and it then
  // contributes to the local.
  std::unique_ptr<Expansion> S_to_M(Point center, SourceRef sources) const;
  std::unique_ptr<Expansion> S_to_L(Point center, SourceRef sources) const;

  //NOTE: These *do* have to wait for the expansion
  // This is a call when on the expansion containted, that will perform the
  // translation, and continue that with the correct code to the expansion LCO
  //
  //TODO: These need more thinking...
  //
  //The results of these are typically added to a particular target or
  // something.
  std::unique_ptr<Expansion> M_to_M(int from_child, double s_size) const;
  std::unique_ptr<Expansion> M_to_L(Index s_index, double s_size,
                                    Index t_index) const;
  std::unique_ptr<Expansion> L_to_L(int to_child, double t_size) const;

  //NOTE: These *do* have to wait
  // This is a call when on the expansion contained, that will perform the
  // translation, and continue that with the correct code to the target LCO
  void M_to_T(TargetRef targets) const;
  void L_to_T(TargetRef targets) const;

  //NOTE: This does not
  //
  void S_to_T(SourceRef sources, TargetRef targets) const;

  //NOTE: This needs to wait on the input expansion
  void add_expansion(const Expansion *temp1);
  //end TODO comment

  std::unique_ptr<Expansion> get_new_expansion(Point center) const;

  //TODO: methods to make the inputs easy
  // these will wrap up the HPX stuff so the user can do "obvious" seeming
  // thigns instead.
  void finalize() const;
  void schedule() const;

 private:
  int type_;
  hpx_addr_t data_;     //this is the LCO
};


ExpansionRef globalize_expansion(std::unique_ptr<Expansion> exp);


} // namespace dashmm


#endif // __DASHMM_EXPANSION_REF_H__
