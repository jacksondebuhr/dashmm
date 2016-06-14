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


/// DistroPolicy is a concept that specifies how the computation represented
/// by the DAG is to be distributed around the available resources. The
/// primary function of this policy is compute_distribution(). These objects
/// can contain parameters that specify the behavior of the policy. These
/// parameters should be specified at construction time.
///
/// DistroPolicy objects should be trivially copyable.


class DistroPolicy {
public:
  /// Computes the distribution of the work represented by the given nodes.
  ///
  /// The only required element of a distribution policy is this one.
  /// This will examine the internal nodes and assign a locality to them.
  /// The sources and targets are assumed to have already had their
  /// locality selected (based on the locality of the data represented by
  /// those nodes of the DAG).
  ///
  /// After this call, all DAGNodes should have their locality set. It is
  /// permitted to parallelize this work, but this routine must not return
  /// before each locality has been set.
  void compute_distribution(DAG &dag);
};
