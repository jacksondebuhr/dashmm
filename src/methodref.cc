/// \file src/methodref.cc
/// \brief Implementation of MethodRef


#include "include/methodref.h"

#include <cassert>

#include <hpx/hpx.h>



namespace dashmm {


/////////////////////////////////////////////////////////////////////
// Interface
/////////////////////////////////////////////////////////////////////


int MethodRef::type() const {
  setup_local_method();
  return met_->type();
}


bool MethodRef::compatible_with(const ExpansionRef expand) const {
  setup_local_method();
  int n_digits = -1; 
  auto locexp = interpret_expansion(expand.type(), nullptr, 0, n_digits);
  return met_->compatible_with(locexp.get());
}


void MethodRef::generate(SourceNode &curr, const ExpansionRef expand) const {
  setup_local_method();
  met_->generate(curr, expand);
}


void MethodRef::aggregate(SourceNode &curr, const ExpansionRef expand) const {
  setup_local_method();
  met_->aggregate(curr, expand);
}


void MethodRef::inherit(TargetNode &curr, const ExpansionRef expand,
                        size_t which_child) const {
  setup_local_method();
  met_->inherit(curr, expand, which_child);
}


void MethodRef::process(TargetNode &curr, std::vector<SourceNode> &consider,
                        bool curr_is_leaf) const {
  setup_local_method();
  met_->process(curr, consider, curr_is_leaf);
}


bool MethodRef::refine_test(bool same_sources_and_targets,
                      const TargetNode &curr,
                      const std::vector<SourceNode> &consider) const {
  setup_local_method();
  return met_->refine_test(same_sources_and_targets, curr, consider);
}


MethodSerial *MethodRef::pin() const {
  if (data_ == HPX_NULL) {
    return nullptr;
  }
  MethodSerial *retval{nullptr};
  hpx_gas_try_pin(data_, (void **)&retval);
  return retval;
}


void MethodRef::unpin() const {
  if (data_ == HPX_NULL) {
    return;
  }
  hpx_gas_unpin(data_);
}


void MethodRef::setup_local_method() const {
  if (met_ != nullptr) {
    return;
  }
  MethodSerial *local = pin();
  met_ = create_method(local->type, local);
  unpin();
}


MethodRef globalize_method(Method *met, hpx_addr_t where) {
  if (met == nullptr) {
    return MethodRef{HPX_NULL};
  }
  MethodSerial *serial = met->release();
  size_t size = serial->size + sizeof(MethodSerial);
  hpx_addr_t data = hpx_gas_alloc_local_at_sync(1, size, 0, where);
  assert(data != HPX_NULL);
  //NOTE: This will be slow as the memory likely did not get allocated as
  // registered memory. However, this will not happen often, so this is perhaps
  // okay.
  hpx_gas_memput_rsync(data, serial, size);
  free(serial);
  return MethodRef{data};
}



} // namespace dashmm
