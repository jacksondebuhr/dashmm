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


#ifndef __DASHMM_METHOD_REF_H__
#define __DASHMM_METHOD_REF_H__


/// \file include/methodref.h
/// \brief Reference object for Methods


#include <vector>

#include <hpx/hpx.h>

#include "include/expansionref.h"
#include "include/method.h"
#include "include/node.h"


namespace dashmm {


/// Reference to a Method
///
/// This object is a reference to an object in GAS. As such, it can be passed
/// by value without worry. The point of the object is to provide an interface
/// to the GAS data without having to worry about the HPX-5 side of things.
///
/// As a reference, it attempts to have as many of the same methods as the
/// underlying Method object.
class MethodRef {
 public:
  explicit MethodRef(hpx_addr_t addr) : data_{addr}, met_{nullptr} { }

  ~MethodRef() {
    if (met_) {
      delete met_;
      met_ = nullptr;
    }
  }

  MethodRef(const MethodRef &other) {
    data_ = other.data();
    met_ = nullptr;
  }

  MethodRef &operator=(MethodRef other) {
    data_ = other.data();
    met_ = nullptr;
    return *this;
  }

  void destroy() {
    if (data_ != HPX_NULL) {
      hpx_gas_free_sync(data_);
    }
  }

  /// The type identifier of the method
  int type() const;

  /// The global address of the stored data
  hpx_addr_t data() const {return data_;}

  /// Determine if this method is compatible with a given expansion
  ///
  /// This will determine the compatability of the method with the given
  /// expansion. This is generally down to deciding if the method uses
  /// local or exponential operators, and checking that the expansion
  /// can provide those operators.
  ///
  /// \param expand - a reference to the expansion in question
  ///
  /// \returns - true in case of compatability; false otherwise
  bool compatible_with(const ExpansionRef expand) const;

  /// Generate the expansion at the leaf of the source tree
  ///
  /// This operation is invoked at the leaves of the source tree to generate
  /// the expansions for the sources. Typically, these will be multipole
  /// expansions.
  ///
  /// \param curr - the current node of the source tree (will be a leaf)
  /// \param expand - a reference to a prototype expansion that can be used
  ///                 to generate the expansion for the given node.
  void generate(SourceNode &curr, const ExpansionRef expand) const;

  /// Combine expansions from children of an internal source node
  ///
  /// This operation is invoked on internal nodes of the source tree to
  /// combine the expansions of the children on the given node into the
  /// expansion for the internal node.
  ///
  /// \param curr - the current node of the source tree (will be internal)
  /// \param expand - a reference to a prototype expansion that can be used
  ///                 to aggregate the expansions of the given node's children.
  void aggregate(SourceNode &curr, const ExpansionRef expand) const;

  /// Inherit an expansion from a target node's parent
  ///
  /// This operation is invoked on nodes in the target tree to inherit the
  /// effect of the expansion collected at the parent of the given node.
  ///
  /// \param curr - the current node of the target tree
  /// \param expand - a reference to a prototype expansion that might be used
  /// \which_child - which child @p curr is of its parent
  void inherit(TargetNode &curr, const ExpansionRef expand,
               size_t which_child) const;

  /// Process the list of source nodes for a given target node
  ///
  /// The bulk of the work in the traversal of the target tree is in process.
  /// This operation is called on each target node and it will look at the
  /// set of source nodes in @p consider, and take the appropriate action
  /// given the usability of those nodes.
  ///
  /// \param curr - the target node in question
  /// \param consider - a vector of the source nodes under consideration
  /// \param curr_is_leaf - indicates if @p curr is a leaf node
  void process(TargetNode &curr, std::vector<SourceNode> &consider,
               bool curr_is_leaf) const;

  /// Decide if the node in the target tree should be refined
  ///
  /// In some advanced methods, it can be desirable to halt refinement of the
  /// target tree before the specified refinement limit is reached. This
  /// operation performs that test.
  ///
  /// \param same_sources_and_targets - are the source and target points the
  ///                                   same
  /// \param curr - the target node being examined
  /// \param consider - the list of source nodes being considered at this node
  ///
  /// \returns - true if the refinement should proceed; false otherwise
  bool refine_test(bool same_sources_and_targets,
                        const TargetNode &curr,
                        const std::vector<SourceNode> &consider) const;

 private:
  MethodSerial *pin() const;
  void unpin() const;
  void setup_local_method() const;

  hpx_addr_t data_;
  mutable Method *met_;
};


/// Create a global version of a given method
///
/// This will export the local data represented by the given method into the
/// global address space, and return a reference to the globalized data. This
/// relies on release() from the specific Method.
///
/// \param met - the method to globalize
/// \param where - an HPX address indicating an address which should be local
///                to the globalized method
///
/// \returns - a reference to the resulting method
MethodRef globalize_method(Method *met, hpx_addr_t where);


} // namespace dashmm


#endif // __DASHMM_METHOD_REF_H__
