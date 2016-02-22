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


/// To qualify for the Method concept in DASHMM, a class must satisfy the
/// following criteria. This file is not included anywhere in DASHMM, but it
/// is in the source distribution as an example, and to explain the
/// Method concept.


/// Methods in DASHMM are template classes parameterized over the types of
/// sources and tagets as well as the Expansion. A full description of the
/// requirements of the Source, Target and Expansion can be found elsewhere.
///
/// When creating a user-defined Method, the name Method in the following
/// should be replaced by the name of the new Method type.
template <typename Source, typename Target,
          template <typename, typename> class Expansion>
class Method {
 public:
  /// These are all useful aliases to define, given the heavily templated
  /// nature of DASHMM.
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Method<Source, Target, Expansion>;
  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, Method>;
  using sourcenode_t = SourceNode<Source, Target, Expansion, Method>;
  using targetnode_t = TargetNode<Source, Target, Expansion, Method>;

  /// TODO: The method will need a default constructor.

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
  /// \param n_digits - the accuracy parameter for the expansion in question.
  void generate(sourcenode_t &curr, int n_digits) const;

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
  /// \param n_digits - the accuracy parameter for the expansion in question.
  void aggregate(sourcenode_t &curr, int n_digits) const;

  /// Inherit an expansion from a target node's parent
  ///
  /// This operation is invoked on nodes in the target tree to inherit the
  /// effect of the expansion collected at the parent of the given node.
  ///
  /// \param curr - the current node of the target tree
  /// \param n_digits - the accuracy parameter for the expansion in question.
  /// \which_child - which child @p curr is of its parent
  void inherit(targetnode_t &curr, int n_digits, size_t which_child) const;

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
  void process(targetnode_t &curr, std::vector<sourcenode_t> &consider,
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
                   const targetnode_t &curr,
                   const std::vector<sourcenode_t> &consider) const;

 private:
  // Any data should be trivially copyable
};
