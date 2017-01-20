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
/// \brief Implementation of JSON format DAG output
///
/// The intent of this file it to make an easily digestible form of the
/// DAG information for use in visualization tools. This implements JSON
/// formatted data, as it is human readable (if boring) and because there
/// are quality JSON readers in a wide array of languages and frameworks.
///
/// NOTE: This is not as robust as other portions of DASHMM, and should be
/// considered to be experimental. The default mode of operation of DASHMM
/// will not even call into these routines, and must be manually enabled
/// by modifying the source code of the library. Please note the intentional
/// vagueness.


#include "dashmm/dag.h"

#include <cstdio>

#include <algorithm>
#include <limits>
#include <map>
#include <string>
#include <utility>

#include "dashmm/index.h"


namespace dashmm {


namespace {


enum class NodeType {
  Source,
  Multipole,
  Target,
  Local,
  Intermediate,
  Unknown
};


class Edge {
 public:
  int source;
  int target;
  int weight;
  Operation op;
};


class Node {
 public:
  std::string id;
  NodeType type;
  int locality;
  int depth;

  Node() : id{""}, type{NodeType::Unknown}, locality{-1}, depth{0} { }
};


std::string node_type_to_print(NodeType type) {
  switch (type) {
  case NodeType::Source:
    return std::string("source");
    break;
  case NodeType::Multipole:
    return std::string("multipole");
    break;
  case NodeType::Target:
    return std::string("target");
    break;
  case NodeType::Local:
    return std::string("local");
    break;
  case NodeType::Intermediate:
    return std::string("Intermediate");
    break;
  case NodeType::Unknown:
    return std::string("ERROR");
  }
  return std::string("ERROR");
}


std::string edge_code_to_print(Operation op) {
  switch (op) {
  case Operation::Nop:
    return std::string("Nop");
    break;
  case Operation::StoM:
    return std::string("StoM");
    break;
  case Operation::StoL:
    return std::string("StoL");
    break;
  case Operation::MtoM:
    return std::string("MtoM");
    break;
  case Operation::MtoL:
    return std::string("MtoL");
    break;
  case Operation::LtoL:
    return std::string("LtoL");
    break;
  case Operation::MtoT:
    return std::string("MtoT");
    break;
  case Operation::LtoT:
    return std::string("LtoT");
    break;
  case Operation::StoT:
    return std::string("StoT");
    break;
  case Operation::MtoI:
    return std::string("MtoI");
    break;
  case Operation::ItoI:
    return std::string("ItoI");
    break;
  case Operation::ItoL:
    return std::string("ItoL");
    break;
  }
  return std::string("ERROR");
}


bool is_edge_from_intermediate(Operation op) {
  return (op == Operation::ItoI || op == Operation::ItoL);
}


void add_dagnode_to_index_entries(std::vector<DAGNode *> &nodes,
                                  std::map<const DAGNode *, int> &dtoi) {
  for (size_t i = 0; i < nodes.size(); ++i) {
    dtoi[nodes[i]] = dtoi[nullptr]++;
  }
}


std::map<const DAGNode *, int> create_dagnode_to_index(DAG &dag) {
  std::map<const DAGNode *, int> retval{};
  retval[nullptr] = 0;    // we use nullptr for the count

  add_dagnode_to_index_entries(dag.source_leaves, retval);
  add_dagnode_to_index_entries(dag.source_nodes, retval);
  add_dagnode_to_index_entries(dag.target_nodes, retval);
  add_dagnode_to_index_entries(dag.target_leaves, retval);

  return retval;
}


bool skip_SandT_operations(Operation op) {
  bool skip{false};
  if (op == Operation::StoT || op == Operation::StoM || op == Operation::StoL
      || op == Operation::LtoT || op == Operation::MtoT) {
    skip = true;
  }
  return skip;
}


void append_out_edges(std::map<const DAGNode *, int> &dtoi,
                      const std::vector<DAGNode *> &nodes,
                      std::vector<Edge> &edges) {
  // loop over the nodes
  for (size_t i = 0; i < nodes.size(); ++i) {
    // loop over the out edges
    std::vector<DAGEdge> &out = nodes[i]->out_edges;
    for (size_t j = 0; j < out.size(); ++j) {
      if (skip_SandT_operations(out[j].op)) continue;
      edges.emplace_back(
        Edge{dtoi[out[j].source], dtoi[out[j].target], out[j].weight, out[j].op}
      );
    }
  }
}


std::vector<Edge> create_edges(std::map<const DAGNode *, int> &dtoi,
                               DAG &dag) {
  std::vector<Edge> retval{};

  append_out_edges(dtoi, dag.source_leaves, retval);
  append_out_edges(dtoi, dag.source_nodes, retval);
  append_out_edges(dtoi, dag.target_nodes, retval);
  // No target_leaves, as they have no out edges

  return retval;
}


std::vector<Node> create_nodes(std::map<const DAGNode *, int> &dtoi,
                               DAG &dag) {
  std::vector<Node> retval(dtoi.size() - 1);

  for (size_t i = 0; i < dag.source_leaves.size(); ++i) {
    retval[dtoi[dag.source_leaves[i]]].type = NodeType::Source;
    retval[dtoi[dag.source_leaves[i]]].locality
        = dag.source_leaves[i]->locality;
  }
  for (size_t i = 0; i < dag.source_nodes.size(); ++i) {
    if (dag.source_nodes[i]->out_edges.size() == 0) {
      retval[dtoi[dag.source_nodes[i]]].type = NodeType::Multipole;
    } else if (is_edge_from_intermediate(dag.source_nodes[i]->out_edges[0].op)) {
      retval[dtoi[dag.source_nodes[i]]].type = NodeType::Intermediate;
    } else {
      retval[dtoi[dag.source_nodes[i]]].type = NodeType::Multipole;
    }
    retval[dtoi[dag.source_nodes[i]]].locality = dag.source_nodes[i]->locality;
  }
  for (size_t i = 0; i < dag.target_nodes.size(); ++i) {
    if (is_edge_from_intermediate(dag.target_nodes[i]->out_edges[0].op)) {
      retval[dtoi[dag.target_nodes[i]]].type = NodeType::Intermediate;
    } else {
      retval[dtoi[dag.target_nodes[i]]].type = NodeType::Local;
    }
    retval[dtoi[dag.target_nodes[i]]].locality = dag.target_nodes[i]->locality;
  }
  for (size_t i = 0; i < dag.target_leaves.size(); ++i) {
    retval[dtoi[dag.target_leaves[i]]].type = NodeType::Target;
    retval[dtoi[dag.target_leaves[i]]].locality
        = dag.target_leaves[i]->locality;
  }
  return retval;
}


void compute_depths(std::vector<Node> &nodes, int &n_src, int &n_trg) {
  int source_num{1};
  int target_num{1};
  for (size_t i = 0; i < nodes.size(); ++i) {
    if (nodes[i].type == NodeType::Source) {
      nodes[i].depth = -source_num;
      ++source_num;
    } else if (nodes[i].type == NodeType::Target) {
      nodes[i].depth = target_num;
      ++target_num;
    }
  }
  n_src = source_num - 1;
  n_trg = target_num - 1;
}


void create_unique_node_ids(std::vector<Node> &nodes) {
  int n_source{0};
  int n_multi{0};
  int n_local{0};
  int n_target{0};
  int n_interm{0};

  for (size_t i = 0; i < nodes.size(); ++i) {
    switch (nodes[i].type) {
    case NodeType::Source:
      nodes[i].id = std::string("S") + std::to_string(n_source++);
      break;
    case NodeType::Multipole:
      nodes[i].id = std::string("M") + std::to_string(n_multi++);
      break;
    case NodeType::Local:
      nodes[i].id = std::string("L") + std::to_string(n_local++);
      break;
    case NodeType::Target:
      nodes[i].id = std::string("T") + std::to_string(n_target++);
      break;
    case NodeType::Intermediate:
      nodes[i].id = std::string("I") + std::to_string(n_interm++);
      break;
    case NodeType::Unknown:
      assert(0);
      break;
    }
  }
}


void print_json_nodes(FILE *ofd, std::vector<Node> &nodes) {
  fprintf(ofd, "  \"nodes\": [\n");

  for (size_t i = 0; i < nodes.size(); ++i) {
    if (i) {
      fprintf(ofd, ",\n");
    }
    fprintf(ofd, "    {\n      \"id\": \"%s\",\n", nodes[i].id.c_str());
    std::string temp = node_type_to_print(nodes[i].type);
    fprintf(ofd, "      \"type\": \"%s\",\n", temp.c_str());
    fprintf(ofd, "      \"locality\": %d,\n      \"depth\": %d\n    }",
            nodes[i].locality, nodes[i].depth);
  }

  fprintf(ofd, "\n  ],\n");
}


void print_json_links(FILE *ofd, std::vector<Edge> &edges) {
  fprintf(ofd, "  \"links\": [\n");

  for (size_t i = 0; i < edges.size(); ++i) {
    if (i) {
      fprintf(ofd, ",\n");
    }
    std::string oper = edge_code_to_print(edges[i].op);
    fprintf(ofd, "    {\"source\": %d, \"target\": %d, \"operation\": \"%s\", \"weight\": %d}",
            edges[i].source, edges[i].target, oper.c_str(), edges[i].weight);
  }

  fprintf(ofd, "\n  ]\n");
}


void print_json(std::vector<Node> &nodes, std::vector<Edge> &edges,
                int ns, int nt, const std::string &fname) {
  FILE *ofd = fopen(fname.c_str(), "w");
  assert(ofd);

  fprintf(ofd, "{\n");
  fprintf(ofd, "  \"nsources\": %d,\n", ns);
  fprintf(ofd, "  \"ntargets\": %d,\n", nt);

  print_json_nodes(ofd, nodes);
  print_json_links(ofd, edges);

  fprintf(ofd, "}");

  fclose(ofd);
}


struct CSVEdge {
  Index s_idx;
  int s_loc;
  Operation op;
  Index t_idx;
  int t_loc;
  DAGNode *source;
  DAGNode *target;
};


void add_out_edges_from_node(DAGNode *node, std::vector<CSVEdge> &edges) {
  for (size_t i = 0; i < node->out_edges.size(); ++i) {
    edges.emplace_back(CSVEdge{node->idx, node->locality,
                               node->out_edges[i].op,
                               node->out_edges[i].target->idx,
                               node->out_edges[i].target->locality,
                               node, node->out_edges[i].target});
  }
}


std::vector<CSVEdge> create_csv_edges(DAG &dag) {
  std::vector<CSVEdge> retval{};

  for (size_t i = 0; i < dag.source_leaves.size(); ++i) {
    add_out_edges_from_node(dag.source_leaves[i], retval);
  }

  for (size_t i = 0; i < dag.source_nodes.size(); ++i) {
    add_out_edges_from_node(dag.source_nodes[i], retval);
  }

  for (size_t i = 0; i < dag.target_nodes.size(); ++i) {
    add_out_edges_from_node(dag.target_nodes[i], retval);
  }

  // No need for targets, as they have no out edges

  return retval;
}


void print_edges_to_file(std::vector<CSVEdge> &edges, std::string fname) {
  FILE *ofd = fopen(fname.c_str(), "w");
  assert(ofd);

  for (size_t i = 0; i < edges.size(); ++i) {
    std::string opstr = edge_code_to_print(edges[i].op);
    fprintf(ofd, "%d %d %d %d - %d - %s - %d %d %d %d - %d - %p %p\n",
            edges[i].s_idx.x(), edges[i].s_idx.y(), edges[i].s_idx.z(),
            edges[i].s_idx.level(), edges[i].s_loc, opstr.c_str(),
            edges[i].t_idx.x(), edges[i].t_idx.y(), edges[i].t_idx.z(),
            edges[i].t_idx.level(), edges[i].t_loc,
            edges[i].source, edges[i].target);
  }

  fclose(ofd);
}


bool csvedge_comp(CSVEdge a, CSVEdge b) {
  if (a.s_idx.level() < b.s_idx.level()) return true;
  if (a.s_idx.x() < b.s_idx.x()) return true;
  if (a.s_idx.y() < b.s_idx.y()) return true;
  if (a.s_idx.z() < b.s_idx.z()) return true;

  if (a.t_idx.level() < b.t_idx.level()) return true;
  if (a.t_idx.x() < b.t_idx.x()) return true;
  if (a.t_idx.y() < b.t_idx.y()) return true;
  if (a.t_idx.z() < b.t_idx.z()) return true;

  if (a.op < b.op) return true;

  return false;
}


} // unnamed namespace


void DAG::toJSON(std::string fname) {
  // Create a mapping from DAGNode * to index
  auto dagnode_to_index = create_dagnode_to_index(*this);

  // Create the edge records
  std::vector<Edge> E = create_edges(dagnode_to_index, *this);

  // Get the nodes made
  std::vector<Node> N = create_nodes(dagnode_to_index, *this);

  int n_s, n_t;
  compute_depths(N, n_s, n_t);
  create_unique_node_ids(N);

  // output
  print_json(N, E, n_s, n_t, fname);
}


void DAG::toEdgeCSV(std::string fname) {
  std::vector<CSVEdge> edges = create_csv_edges(*this);
  std::stable_sort(edges.begin(), edges.end(), csvedge_comp);
  print_edges_to_file(edges, fname);
}


size_t DAG::node_count() const {
  return (source_leaves.capacity() + source_nodes.capacity()
          + target_nodes.capacity() + target_leaves.capacity());
}


size_t DAG::edge_count() const {
  size_t retval{0};

  for (auto i = source_leaves.begin(), e = source_leaves.end(); i != e; ++i) {
    retval += (*i)->out_edges.capacity() + (*i)->in_edges.capacity();
  }
  for (auto i = source_nodes.begin(), e = source_nodes.end(); i != e; ++i) {
    retval += (*i)->out_edges.capacity() + (*i)->in_edges.capacity();
  }
  for (auto i = target_nodes.begin(), e = target_nodes.end(); i != e; ++i) {
    retval += (*i)->out_edges.capacity() + (*i)->in_edges.capacity();
  }
  for (auto i = target_leaves.begin(), e = target_leaves.end(); i != e; ++i) {
    retval += (*i)->out_edges.capacity() + (*i)->in_edges.capacity();
  }

  return retval;
}


} // dashmm
