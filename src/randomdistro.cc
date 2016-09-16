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


#include "builtins/randomdistro.h"

#include <cassert>

#include <random>
#include <vector>


namespace dashmm {


void RandomDistro::compute_distribution(DAG &dag) {
  // check that all sources and targets have locality set
  for (size_t i = 0; i < dag.source_leaves.size(); ++i) {
    assert(dag.source_leaves[i]->locality != -1);
  }
  for (size_t i = 0; i < dag.target_leaves.size(); ++i) {
    assert(dag.target_leaves[i]->locality != -1);
  }

  // set up RNG
  std::mt19937 engine(seed_);
  std::uniform_int_distribution<int> uniform(0, hpx_get_num_ranks() - 1);

  // loop over the internals and set at random
  for (size_t i = 0; i < dag.source_nodes.size(); ++i) {
    if (dag.source_nodes[i]->locality != -1) continue;

    dag.source_nodes[i]->locality = uniform(engine);
  }
  for (size_t i = 0; i < dag.target_nodes.size(); ++i) {
    if (dag.target_nodes[i]->locality != -1) continue;

    dag.target_nodes[i]->locality = uniform(engine);
  }
}


} // dashmm
