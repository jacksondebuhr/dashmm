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


/// \file
/// \brief Implement DAG objects


#include "dashmm/dag.h"

#include <cstdio>

#include <vector>


namespace dashmm {

namespace {
  using edge_table_t = std::vector<std::vector<std::vector<int>>>;

  int optoint(Operation op) {
    switch(op) {
      case Operation::Nop:
        return 0;
        break;
      case Operation::StoM:
        return 1;
        break;
      case Operation::StoL:
        return 2;
        break;
      case Operation::MtoM:
        return 3;
        break;
      case Operation::MtoL:
        return 4;
        break;
      case Operation::LtoL:
        return 5;
        break;
      case Operation::MtoT:
        return 6;
        break;
      case Operation::LtoT:
        return 7;
        break;
      case Operation::StoT:
        return 8;
        break;
      case Operation::MtoI:
        return 9;
        break;
      case Operation::ItoI:
        return 10;
        break;
      case Operation::ItoL:
        return 11;
        break;
    }
    return 0;
  }

  edge_table_t zero_table(int n) {
    edge_table_t edges(
        12, std::vector<std::vector<int>>(
            n, std::vector<int>(n, 0)));
    return edges;
  }

  void collect_edges(edge_table_t &edges, const std::vector<DAGNode *> &nodes) {
    for (size_t i = 0; i < nodes.size(); ++i) {
      for (size_t j = 0; j < nodes[i]->out_edges.size(); ++j) {
        DAGEdge &e = nodes[i]->out_edges[j];
        edges[optoint(e.op)][nodes[i]->locality][e.target->locality] += 1;
      }
    }
  }
}


Index DAGNode::index() const {
  return parent_->index();
}

bool DAGNode::is_parts() const {
  return parent_->parts() == this;
}

bool DAGNode::is_normal() const {
  return parent_->normal() == this;
}

bool DAGNode::is_interm() const {
  return parent_->interm() == this;
}

void *DAGNode::tree_node() const {
  return parent_->tree_node();
}

void DAG::printedges(int n) {
  edge_table_t edges = zero_table(n);
  collect_edges(edges, source_leaves);
  collect_edges(edges, source_nodes);
  collect_edges(edges, target_nodes);

  FILE *ofd = fopen("dag.txt", "w");
  for (int i = 0; i < 12; ++i) {
    fprintf(ofd, "Table for operation %d:\n", i);
    fprintf(ofd, "------------------------------\n");
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        fprintf(ofd, "%d -> %d : %d\n", j, k, edges[i][j][k]);
      }
    }
    fprintf(ofd, "\n");
  }
  fclose(ofd);
}


} // dashmm


