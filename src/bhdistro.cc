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


#include "builtins/bhdistro.h"

#include <vector>


namespace dashmm {

void BHDistro::compute_distribution(DAG &dag) {
  std::queue<DAGNode *> nodes = collect_readies(dag);

  while (!nodes.empty()) {
    DAGNode *curr = nodes.front();
    compute_locality(curr);
    mark_upstream_nodes(curr, nodes);
    nodes.pop();
  }
}


// Collect all nodes of the DAG which have no out edges. These are either
// sources, which already have a locality, or are nodes that do not matter
// much, in which case a random assignment is okay.
std::queue<DAGNode *> BHDistro::collect_readies(DAG &dag) {
  std::queue<DAGNode *> retval{};

  for (size_t i = 0; i < dag.source_nodes.size(); ++i) {
    if (dag.source_nodes[i]->out_edges.size() == 0) {
      retval.push(dag.source_nodes[i]);
    }
  }

  for (size_t i = 0; i < dag.target_nodes.size(); ++i) {
    if (dag.target_nodes[i]->out_edges.size() == 0) {
      retval.push(dag.target_nodes[i]);
    }
  }

  for (size_t i = 0; i < dag.target_leaves.size(); ++i) {
    assert(dag.target_leaves[i]->out_edges.size() == 0);
    retval.push(dag.target_leaves[i]);
  }

  return retval;
}


void BHDistro::compute_locality(DAGNode *node) {
  // It already has a locality
  if (node->locality >= 0) return;

  // A rare case where a node does not enter into the computation. So pick
  // at random. A better idea might be to wait until all others are placed
  // ignorning this node, and then pick the best given the localities of the
  // upstream nodes, but for now, we do something simple.
  if (node->out_edges.size() == 0) {
    node->locality = 0;
  }

  // The typical case; count up weights to each locality
  int n_ranks = hpx_get_num_ranks();
  std::vector<int> bins(n_ranks, 0);
  for (size_t i = 0; i < node->out_edges.size(); ++i) {
    int loc = node->out_edges[i].target->locality;
    assert(loc >= 0 && loc < n_ranks);
    bins[loc] += node->out_edges[i].weight;
  }

  // Find max - currently, the lowest locality in a tie will win.
  int maxval{bins[0]};
  size_t maxidx{0};
  for (int i = 1; i < n_ranks; ++i) {
    if (bins[i] > maxval) {
      maxval = bins[i];
      maxidx = i;
    }
  }

  // Set this node to have locality that matches the bulk of outgoing work
  node->locality = maxidx;
}


void BHDistro::mark_upstream_nodes(DAGNode *node,
                                   std::queue<DAGNode *> &master) {
  assert(node->locality >= 0);  // Need to be sure that the locality is set.

  for (size_t i = 0; i < node->in_edges.size(); ++i) {
    node->in_edges[i].source->color += 1;
    if ((size_t)node->in_edges[i].source->color
            == node->in_edges[i].source->out_edges.size()) {
      master.push(node->in_edges[i].source);
    }
  }
}


bool BHDistro::distribution_complete(DAG &dag) {
  for (size_t i = 0; i < dag.source_leaves.size(); ++i) {
    if (dag.source_leaves[i]->locality < 0) {
      return false;
    }
  }

  for (size_t i = 0; i < dag.source_nodes.size(); ++i) {
    if (dag.source_nodes[i]->locality < 0) {
      return false;
    }
  }

  for (size_t i = 0; i < dag.target_nodes.size(); ++i) {
    if (dag.target_nodes[i]->locality < 0) {
      return false;
    }
  }

  for (size_t i = 0; i < dag.target_leaves.size(); ++i) {
    if (dag.target_leaves[i]->locality < 0) {
      return false;
    }
  }

  return true;
}


} // namespace dashmm
