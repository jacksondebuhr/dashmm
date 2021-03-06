// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2017, Trustees of Indiana University,
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
///
/// Methods in DASHMM are template classes parameterized over the types of
/// sources and targets, and the Expansion. A full description of the
/// requirements of the Source, Target and Expansion can be found elsewhere.
///
/// When creating a user-defined Method, the name Method in the following
/// should be replaced by the name of the new Method type.
///
/// Methods operate on the DAG which is explicitly represented in the DAGInfo
/// and DAGNode objects. See their documentation for more details.
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
  using sourcenode_t = Node<Source>;
  using targetnode_t = Node<Target>;

  /// Each method requires the specification of a distribution policy aliased
  /// to distropolicy_t. The best all around distribution policy included with
  /// DASHMM will be made available through dashmm/defaultpolicy.h as
  /// DefaultDistributionPolicy. If the implementer of a method has a better
  /// notion about how to distribute the nodes of the DAG for the particlar
  /// method, then a policy can be defined, and set as the new method's
  /// policy here.
  using distropolicy_t = DefaultDistributionPolicy;

  /// Methods require a default constructor.
  Method();

  /// Generate the expansion at the leaf of the source tree
  ///
  /// This operation is invoked at the leaves of the source tree to generate
  /// links in the DAG from sources into internal nodes of the DAG. Typically
  /// this will schedule S->M operations.
  ///
  /// \param curr - the current node of the source tree (will be a leaf)
  /// \param domain - the domain geometry for the tree
  void generate(sourcenode_t *curr, DomainGeometry *domain) const;

  /// Combine expansions from children of an internal source node
  ///
  /// This operation is invoked on internal nodes of the source tree to
  /// combine add the edges in the DAG that connect multipole moments on the
  /// source tree into other multipole moments.
  ///
  /// \param curr - the current node of the source tree (will be internal)
  /// \param domain - the domain geometry for the tree
  void aggregate(sourcenode_t *curr, DomainGeometry *domain) const;

  /// Inherit an expansion from a target node's parent
  ///
  /// This operation is invoked on nodes in the target tree to inherit the
  /// effect of the expansion collected at the parent of the given node.
  ///
  /// \param curr - the current node of the target tree
  /// \param domain - the domain geometry for the tree
  /// \param curr_us_leaf - indicates the the current node is a leaf
  void inherit(targetnode_t *curr, DomainGeometry *domain,
               bool curr_is_leaf) const;

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
  /// \param domain - the domain geometry for the tree
  void process(targetnode_t *curr, std::vector<sourcenode_t *> &consider,
               bool curr_is_leaf, DomainGeometry *domain) const;

  /// Decide if the node in the target tree should be refined
  ///
  /// In some advanced methods, it can be desirable to halt refinement of the
  /// target tree before the specified refinement limit is reached. This
  /// operation performs that test.
  ///
  /// NOTE: In the distributed case, the nodes may not have information about
  /// the number of sources or targets in the node. So any refinement test
  /// cannot depend on those data. Further, though this is initially
  /// conceptually about the refinement of the tree, in distributed cases,
  /// to allow for consistent views on each rank, the refinement proceeds
  /// anyway. Instead, this is more to do with when the DAG stops being
  /// discovered in the tree.
  ///
  /// \param same_sources_and_targets - are the source and target points the
  ///                                   same
  /// \param curr - the target node being examined
  /// \param consider - the list of source nodes being considered at this node
  ///
  /// \returns - true if the refinement should proceed; false otherwise
  bool refine_test(bool same_sources_and_targets,
                   const targetnode_t *curr,
                   const std::vector<sourcenode_t *> &consider) const;

 private:
  // Any data should be trivially copyable
};
