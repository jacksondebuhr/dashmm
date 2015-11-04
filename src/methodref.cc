#include "include/methodref.h"

#include <cassert>

#include <hpx/hpx.h>

//other dashmm


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
  return met_->compatible_with(expand);
}


void MethodRef::generate(SourceNode &curr, const ExpansionRef expand) const {
  setup_local_method();
  met_->generate(curr, expand);
}


void MethodRef::aggregate(SourceNode &curr, const ExpansionRef expand) const {
  setup_local_method();
  met_->generate(curr, expand);
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


void bool MethodRef::refine_test(bool same_sources_and_targets,
                      const TargetNode &curr,
                      const std::vector<SourceNode> &consider) const {
  setup_local_method();
  return met_->refine_test(same_sources_and_targets, curr, consider);
}


MethodSerial *MethodRef::pin() {
  if (data_ == HPX_NULL) {
    return nullptr;
  }
  MethodSerial *retval{nullptr};
  hpx_gas_try_pin(data_, (void **)&retval);
  return retval;
}


void MethodRef::unpin() {
  if (data_ == HPX_NULL) {
    return;
  }
  hpx_gas_unpin(data_);
}


void MethodRef::setup_local_method() {
  if (met_ != nullptr) {
    return;
  }
  MethodSerial *local = pin();
  met_ = create_method(local->type, local->size,
                      static_cast<void *>local->data);
  unpin();
}


MethodRef globalize_method(Method *met, hpx_addr_t where) {
  if (met == nullptr) {
    return HPX_NULL;
  }
  MethodSerialPtr serial = met->serialize(true);
  size_t size = serial->size + sizeof(MethodSerial);
  hpx_addr_t data = hpx_gas_alloc_local_at_sync(size, 0, where);
  assert(data != HPX_NULL);
  hpx_gas_memput_rsync(data, serial.get(), size);
  return MethodRef{data};
}



} // namespace dashmm
