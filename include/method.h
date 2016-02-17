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


#ifndef __DASHMM_METHOD_H__
#define __DASHMM_METHOD_H__


/// \file include/method.h
/// \brief Abstract interface for Method objects


#include <vector>


#include "include/expansionlco.h"
#include "include/node.h"
#include "include/types.h"


namespace dashmm {


/// Abstract interface for Methods used in DASHMM
///
/// This interface specifies the requirements for methods that a user
/// may add to DASHMM. In general, the user will need to implement the method
/// and then somewhere in their program call register_user_method() providing
/// the creation functions that the user also implements. For details, see the
/// DASHMM Advanced User Guide.
template <typename Source, typename Target,
          template <typename, typename> class Expansion>
class Method {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;

  using method_t = Method<Source, Target, Expansion>;

  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, Method>;

  /// Generate the expansion at the leaf of the source tree
  ///
  /// This operation is invoked at the leaves of the source tree to generate
  /// the expansions for the sources. Typically, these will be multipole
  /// expansions.
  ///
  /// It is assumed that generate will either create the expansion with the
  /// needed data, or it will create an empty expansion and schedule the
  /// contribution to the expansion. Internally, DASHMM will call finalize on
  /// the expansion after generate() is called for a particular node.
  /// Further, generate will not return until the expansion for @p curr has
  /// been set. This does not require that all contributions have been made to
  /// that expansion.
  ///
  /// \param curr - the current node of the source tree (will be a leaf)
  /// \param expand - a reference to a prototype expansion that can be used
  ///                 to generate the expansion for the given node.
  void generate(SourceNode &curr, const expansionlco_t expand) const;

  /// Combine expansions from children of an internal source node
  ///
  /// This operation is invoked on internal nodes of the source tree to
  /// combine the expansions of the children on the given node into the
  /// expansion for the internal node.
  ///
  /// It is assumed that aggregate creates the expansion and schedules all
  /// contributions to that expansion.
  ///
  /// \param curr - the current node of the source tree (will be internal)
  /// \param expand - a reference to a prototype expansion that can be used
  ///                 to aggregate the expansions of the given node's children.
  void aggregate(SourceNode &curr, const expansionlco_t expand) const;

  /// Inherit an expansion from a target node's parent
  ///
  /// This operation is invoked on nodes in the target tree to inherit the
  /// effect of the expansion collected at the parent of the given node.
  ///
  /// \param curr - the current node of the target tree
  /// \param expand - a reference to a prototype expansion that might be used
  /// \which_child - which child @p curr is of its parent
  void inherit(TargetNode &curr, const expansionlco_t expand,
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
  // Any data should be trivially copyable
};


} // namespace dashmm


#endif // __DASHMM_METHOD_H__
