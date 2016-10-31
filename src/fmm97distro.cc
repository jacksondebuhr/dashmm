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

#include "builtins/fmm97distro.h"

#include <algorithm>
#include <map>
#include <limits>


namespace dashmm {


void FMM97Distro::assign_for_source(DAGInfo &dag, int locality, int height) {
  DAGNode *normal = dag.normal();
  if (normal != nullptr) {
    normal->locality = locality;
    normal->color = height + 1;
  }

  DAGNode *interm = dag.interm();
  if (interm != nullptr) {
    interm->locality = locality;
    interm->color = height + 2;
  }
}

void FMM97Distro::assign_for_target(DAGInfo &dag, int locality) {
  DAGNode *normal = dag.normal();
  DAGNode *interm = dag.interm();

  if (normal != nullptr) {
    normal->locality = locality;
  }

  if (interm != nullptr) {
    // Categorize incoming edges
    std::map<int, int> color;
    std::map<int, int> weight;
    int in_weight = 0;

    for (size_t i = 0; i < interm->in_edges.size(); ++i) {
      int w = interm->in_edges[i].weight;
      int c = interm->in_edges[i].source->color;
      int source_locality = interm->in_edges[i].source->locality;

      color[source_locality] = std::max(color[source_locality], c);
      weight[source_locality] += w;
      in_weight += w;
    }

    int min_weight = std::numeric_limits<int>::max();
    int max_color = std::numeric_limits<int>::min();
    int interm_locality = -1;

    for (auto i = weight.begin(); i != weight.end(); ++i) {
      int source_locality = i->first;
      int w = in_weight - i->second;
      int c = color[source_locality];

      if (w < min_weight) {
        min_weight = w;
        max_color = c;
        interm_locality = source_locality;
      } else if (w == min_weight && c > max_color) {
        max_color = c;
        interm_locality = source_locality;
      }
    }

    interm->locality = interm_locality;
    assert(interm_locality != -1);
  }
}


void FMM97Distro::compute_distribution(DAG &dag) {
  // TODO:
  // There may exist some nodes in the source_nodes container that do not have
  // outgoing edges. In the future, these nodes can be eliminated (no allocation
  // or work scheduled).
}


} // dashmm
