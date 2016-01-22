// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license. See the LICENSE file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================


#ifndef __DASHMM_DIRECT_METHOD_H__
#define __DASHMM_DIRECT_METHOD_H__


/// \file include/direct_method.h
/// \brief Direction summation method


#include "include/ids.h"
#include "include/expansionref.h"
#include "include/method.h"
#include "include/node.h"


namespace dashmm {

//
class Direct : public Method {
 public:
  Direct() : local_{nullptr} {
    local_ = static_cast<MethodSerial *>(malloc(sizeof(MethodSerial)));
    assert(local_);
    local_->type = kMethodDirect;
    local_->size = 0;
  }

  ~Direct() {
    if (local_) {
      free(local_);
      local_ = nullptr;
    }
  }

  MethodSerial *release() override {
    MethodSerial *retval = local_;
    local_ = nullptr;
    return retval;
  }

  size_t bytes() const override {
    return sizeof(MethodSerial);
  }

  int type() const override {return kMethodDirect;}

  bool compatible_with(const Expansion *expand) const override {return true;}

  void generate(SourceNode &curr, const ExpansionRef expand) const override;
  void aggregate(SourceNode &curr, const ExpansionRef expand) const override;
  void inherit(TargetNode &curr, const ExpansionRef expand,
               size_t which_child) const override;
  void process(TargetNode &curr, std::vector<SourceNode> &consider,
               bool curr_is_leaf) const override;

  //One option would be to cause this test to always return false, meaning the
  // trees would not refine at all, and thus create two huge nodes, and then
  // run S->T for the given expansion. However, we can get some parallelism if
  // we go ahead and refine where we can. Just because we are using a dumb
  // method, does not mean we have to do that dumb thing in serial.
  bool refine_test(bool same_sources_and_targets, const TargetNode &curr,
                   const std::vector<SourceNode> &consider) const override {
    return true;
  }

 private:
  MethodSerial *local_;
};


} // namespace dashmm


#endif // __DASHMM_DIRECT_METHOD_H__
