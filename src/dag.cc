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


/// \file src/dag.cc
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
/// by defining the correct symbol at compile time. Note that this is left
/// vague on purpose.


#include "dashmm/dag.h"

#include <cstdio>

#include <map>
#include <string>

#include "dashmm/index.h"


namespace dashmm {


namespace {


enum class NodeType {
  Source,
  Multipole,
  Target,
  Local,
  Unknown
};


class Edge {
 public:
  int source;
  int target;
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


void append_out_edges(std::map<const DAGNode *, int> &dtoi,
                      const std::vector<DAGNode *> &nodes,
                      std::vector<Edge> &edges) {
  // loop over the nodes
  for (size_t i = 0; i < nodes.size(); ++i) {
    // loop over the out edges
    std::vector<DAGEdge> &out = nodes[i]->out_edges;
    for (size_t j = 0; j < out.size(); ++j) {
      edges.emplace_back(
        Edge{dtoi[out[j].source], dtoi[out[j].target], out[j].op}
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
    retval[dtoi[dag.source_nodes[i]]].type = NodeType::Multipole;
    retval[dtoi[dag.source_nodes[i]]].locality = dag.source_nodes[i]->locality;
  }
  for (size_t i = 0; i < dag.target_nodes.size(); ++i) {
    retval[dtoi[dag.target_nodes[i]]].type = NodeType::Local;
    retval[dtoi[dag.target_nodes[i]]].locality = dag.target_nodes[i]->locality;
  }
  for (size_t i = 0; i < dag.target_leaves.size(); ++i) {
    retval[dtoi[dag.target_leaves[i]]].type = NodeType::Target;
    retval[dtoi[dag.target_leaves[i]]].locality
        = dag.target_leaves[i]->locality;
  }
  return retval;
}


// TODO make this actually do the depth
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
    fprintf(ofd, "    {\"source\": %d, \"target\": %d, \"operation\": \"%s\"}",
            edges[i].source, edges[i].target, oper.c_str());
  }

  fprintf(ofd, "\n  ]\n");
}


void print_json(std::vector<Node> &nodes, std::vector<Edge> &edges,
                int ns, int nt) {
  FILE *ofd = fopen("dag.json", "w");
  assert(ofd);

  fprintf(ofd, "{\n");
  fprintf(ofd, "  \"nsources\": %d,\n", ns);
  fprintf(ofd, "  \"ntargets\": %d,\n", nt);

  print_json_nodes(ofd, nodes);
  print_json_links(ofd, edges);

  fprintf(ofd, "}");

  fclose(ofd);
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
  print_json(N, E, n_s, n_t);
}


} // namespace dashmm
