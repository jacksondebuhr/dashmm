// =============================================================================
//  This file is part of:
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  DASHMM is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  DASHMM is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with DASHMM. If not, see <http://www.gnu.org/licenses/>.
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
///
/// Some incoming nodes of the DAG will already have a locality set. Only
/// those nodes that have locality of -1 initially should be given a
/// locality during the distribution computation.
class DistroPolicy {
public:
  /// Construct a default version of the distribution policy. This can
  /// either be a constructor taking no arguments, or a constructor with all
  /// arguments having a default value.
  DistroPolicy();

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
  ///
  /// This routine should be deterministic, and further, the results should
  /// not depend on the precise ordering of DAG nodes in @p dag, or the edges
  /// in the nodes of the DAG.
  void compute_distribution(DAG &dag);
};
