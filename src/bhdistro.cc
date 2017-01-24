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


/// \file src/bhdistro.cc
/// \brief Implement BHDistro


#include "builtins/bhdistro.h"

#include <algorithm>


namespace dashmm {


void BHDistro::compute_distribution(DAG &dag) {
  assign_colors_and_sort(dag);
  assign_localities(dag.target_nodes);
  assign_localities(dag.source_nodes);
}


void BHDistro::assign_colors_and_sort(DAG &dag) {
  for (size_t i = 0; i < dag.source_nodes.size(); ++i) {
    DAGNode *node = dag.source_nodes[i];
    node->color = -2 * node->index().level() + (node->is_interm() ? 1 : 0);
  }
  for (size_t i = 0; i < dag.target_nodes.size(); ++i) {
    DAGNode *node = dag.target_nodes[i];
    node->color = 2 * (node->index().level() + 1)
                  - (node->is_interm() ? 1 : 0);
  }
  std::sort(dag.source_nodes.begin(), dag.source_nodes.end(),
            color_comparison);
  std::sort(dag.target_nodes.begin(), dag.target_nodes.end(),
            color_comparison);
}


bool BHDistro::color_comparison(const DAGNode *a, const DAGNode *b) {
  return a->color > b->color;
}

void BHDistro::assign_localities(std::vector<DAGNode *> &nodes) {
  // Upon entry, the nodes are sorted highest to lowest color
  for (size_t i = 0; i < nodes.size(); ++i) {
    compute_locality(nodes[i]);
  }
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
  std::vector<int> bins(n_ranks, 0);    // Is there a better choice here?
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
