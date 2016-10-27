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
