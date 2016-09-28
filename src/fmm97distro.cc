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

void FMM97Distro::compute_distribution(DAG &dag) {
  // color all the DAG node
  for (size_t i = 0; i < dag.source_leaves.size(); ++i) {
    DAGNode *s = dag.source_leaves[i];
    s->color = 0;
    color(s);
  }

  // Confine M and I of the source tree
  for (size_t i = 0; i < dag.source_leaves.size(); ++i) {
    DAGNode *n = dag.source_leaves[i];
    assert(n->locality != -1);
    confine(n, 's');
  }

  // Confine L of the target tree
  for (size_t i = 0; i < dag.target_leaves.size(); ++i) {
    DAGNode *n = dag.target_leaves[i];
    assert(n->locality != -1);
    confine(n, 't');
  }
  
  // Make decision on I of the target tree
  for (size_t i = 0; i < dag.target_nodes.size(); ++i) {
    DAGNode *n = dag.target_nodes[i];
    if (n->locality == -1) {
      assign(n);
    }
  }
}

void FMM97Distro::color(DAGNode *s) {
  for (size_t i = 0; i < s->out_edges.size(); ++i) {
    DAGNode *t = s->out_edges[i].target;
    t->color = std::max(t->color, s->color + 1);
    color(t);
  }
}

void FMM97Distro::confine(DAGNode *n, char type) {
  assert(type == 's' || type == 't');

  if (type == 's') {
    for (size_t i = 0; i < n->out_edges.size(); ++i) {
      DAGNode *target = n->out_edges[i].target; 
      Operation op = n->out_edges[i].op; 
      bool terminate = true; 

      if (op == Operation::MtoI) {
        target->locality = n->locality; 
      } 

      if (op == Operation::MtoM) {
        if (n->locality > target->locality) {
          target->locality = n->locality; 
          terminate = false;
        }
      }

      if (op == Operation::StoM) 
        terminate = false; 

      if (terminate == false) 
        confine(target, type); 
    }
  } 

  if (type == 't') {
    for (size_t i = 0; i < n->in_edges.size(); ++i) {
      DAGNode *source = n->in_edges[i].source;
      Operation op = n->in_edges[i].op;
      bool terminate = true; 

      if (op == Operation::LtoT) 
        terminate = false; 

      if (op == Operation::LtoL) {
        if (n->locality > source->locality) {
          source->locality = n->locality; 
          terminate = false; 
        }
      }

      if (terminate == false) 
        confine(source, type); 
    }
  }
}


void FMM97Distro::assign(DAGNode *n) {
  /// \param n is the expansion LCO for an intermediate expansion of the target
  /// tree. 

  // Categorize incoming edges 
  std::map<int, int> color; 
  std::map<int, int> weight; 
  int in_weight = 0; 

  for (size_t i = 0; i < n->in_edges.size(); ++i) {
    int w = n->in_edges[i].weight; 
    int c = n->in_edges[i].source->color; 
    int source_locality = n->in_edges[i].source->locality; 

    color[source_locality] = std::max(color[source_locality], c); 
    weight[source_locality] += w; 
    in_weight += w;
  }

  int min_weight = std::numeric_limits<int>::max(); 
  int max_color = std::numeric_limits<int>::min(); 
  int locality = -1; 

  for (auto i = weight.begin(); i != weight.end(); ++i) {
    int source_locality = i->first; 

    // Adjust incoming communication cost if \param n is placed on
    // source_locality 
    int w = in_weight - i->second; 

    // Adjust ItoL cost if \param n is placed on source_locality 
    for (size_t j = 0; j < n->out_edges.size(); ++j) {
      int target_locality = n->out_edges[j].target->locality; 
      if (target_locality != source_locality) 
        w += n->out_edges[0].weight; 
    }

    int c = color[source_locality]; 

    if (w < min_weight) {
      min_weight = w; 
      max_color = c; 
      locality = source_locality;
    } else if (w == min_weight && c > max_color) {
      // Increase slack time 
      max_color = c; 
      locality = source_locality; 
    }
  }

  n->locality = locality; 
  assert(locality != -1); 
}

} // dashmm
