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


// This will lay out the requirements on the DistroPolicy object type.
// New policies will need to conform to this.

class DistroPolicy {
public:
  // Distribution policies might need some input parameters if their behavior
  // is tunable in any way. TODO: work out what that means for the use of this
  // policy. Only the tree is going to use it really. Do we also pass in a
  // default argument to evaluate (the policy object setup how it should be)?

  // TODO does this need anything else?
  //
  // NOTE: this will compute the distribution of the non-terminal nodes.
  // That is, anything not a source or target. It is assumed that those nodes
  // have already had their locality computed.
  //
  // Here the source and target nodes are separated out to make some things
  // simpler elsewhere in the library. But it could be useful to have a handle
  // one starting and ending points for a specific policy.
  void compute_distribution(const SharedData<DomainGeometry> &domain,
                            const std::vector<DAGNode *> &sources,
                            const std::vector<DAGNode *> &targets,
                            const std::vector<DAGNode *> &internal);

  // Is this policy going to interact with the execution anywhere else?
  // Perhaps eventually this can control the placment of the sources and
  // targets?

  // NOTE: the DistroPolicy will need the appropriate copy operations
};
