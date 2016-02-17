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


#ifndef __DASHMM_BH_METHOD_H__
#define __DASHMM_BH_METHOD_H__


/// \file include/bh_method.h
/// \brief Declaration of BHMethod


#include "include/expansionlco.h"
#include "include/method.h"
#include "include/node.h"


namespace dashmm {


/// A Method to implement classic Barnes-Hut
///
/// It uses the simple critical angle criterion to decide if a given expansion
/// is usable. There is little that needs explanation in this class, as the
/// methods are essentially just those of the abstract base class Method.
class BH : public Method {
 public:
  BH(double theta) : local_{nullptr} {
    local_ = reinterpret_cast<MethodSerial *>(
            new char [sizeof(MethodSerial) + sizeof(double)]);
    assert(local_);
    local_->type = kMethodBH;
    local_->size = sizeof(double);
    local_->data[0] = theta;
  }

  ~BH() {
    if (local_) {
      delete [] local_;
      local_ = nullptr;
    }
  }

  MethodSerial *release() override {
    MethodSerial *retval{local_};
    local_ = nullptr;
    return retval;
  }

  size_t bytes() const override {
    return (sizeof(MethodSerial) + sizeof(double));
  }

  int type() const override {return kMethodBH;}

  bool compatible_with(const Expansion *expand) const override {return true;}

  void generate(SourceNode &curr, const ExpansionRef expand) const override;
  void aggregate(SourceNode &curr, const ExpansionRef expand) const override;
  void inherit(TargetNode &curr, const ExpansionRef expand,
               size_t which_child) const override;
  void process(TargetNode &curr, std::vector<SourceNode> &consider,
                       bool curr_is_leaf) const override;

  //BH always calls for refinement
  bool refine_test(bool same_sources_and_targets, const TargetNode &curr,
                   const std::vector<SourceNode> &consider) const override {
    return true;
  }

  /// Decide on the usability of an expansion
  ///
  /// This performs the traditional critical angle comparison.
  ///
  /// \param exp_point - the position of the expansion center
  /// \param size - the size of the containing node
  /// \param pos - the position of the point to check for usability.
  ///
  /// \returns - true if the expansion is usable; false otherwise
  bool MAC(Point exp_point, double size, Point pos) const;

  /// Finds the point nearest the given point in a given node
  ///
  /// This is used to compute the point in a target node that is nearest
  /// the source expansion center inside a given target node. This allows the
  /// MAC to be used once to verify for all possible points inside the target
  /// node.
  ///
  /// \param scenter - the center of the source expansion
  /// \param tcenter - the center of the target node
  /// \param tsize - the size of the target node
  ///
  /// \returns - the point in the target node closest to the expansion center
  Point nearest(Point scenter, Point tcenter, double tsize) const;

 private:
  MethodSerial *local_;
};


} // namespace dashmm


#endif // __DASHMM_BH_METHOD_H__
