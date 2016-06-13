#include "dashmm/dagtojson.h"

#include <cstdio>

#include <map>
#include <string>

#include "dashmm/index.h"

namespace dashmm {

namespace {

class FullDAGEdge {
 public:
  const DAGNode *source;
  const DAGNode *target;
  Operation op;
};

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
  }
  return std::string("ERROR");
}


NodeType edge_type_to_node_type(Operation op) {
  switch (op) {
    case Operation::Nop:
      return NodeType::Unknown;
      break;
    case Operation::StoM:
      return NodeType::Source;
      break;
    case Operation::StoL:
      return NodeType::Source;
      break;
    case Operation::MtoM:
      return NodeType::Multipole;
      break;
    case Operation::MtoL:
      return NodeType::Multipole;
      break;
    case Operation::LtoL:
      return NodeType::Local;
      break;
    case Operation::MtoT:
      return NodeType::Multipole;
      break;
    case Operation::LtoT:
      return NodeType::Local;
      break;
    case Operation::StoT:
      return NodeType::Source;
      break;
  }
  return NodeType::Unknown;
}


void collect_full_edges_from_vector(std::vector<DAGNode *> &nodes,
                                    std::vector<FullDAGEdge> &edges) {
  for (size_t i = 0; i < nodes.size(); ++i) {
    for (size_t j = 0; j < nodes[i]->edges.size(); ++j) {
      edges.emplace_back(
        FullDAGEdge{nodes[i], nodes[i]->edges[j].target, nodes[i]->edges[j].op}
      );
    }
  }
}


std::vector<FullDAGEdge> collect_full_dag_edges(
    std::vector<DAGNode *> &source, std::vector<DAGNode *> &target,
    std::vector<DAGNode *> &internal) {
  std::vector<FullDAGEdge> retval{};

  collect_full_edges_from_vector(source, retval);
  collect_full_edges_from_vector(internal, retval);
  collect_full_edges_from_vector(target, retval);

  return retval;
}


void add_dagnode_to_index_entries(std::vector<DAGNode *> &nodes,
                                  std::map<const DAGNode *, int> &dtoi) {
  for (size_t i = 0; i < nodes.size(); ++i) {
    dtoi[nodes[i]] = dtoi[nullptr]++;
  }
}

std::map<const DAGNode *, int> create_dagnode_to_index(
    std::vector<DAGNode *> &source, std::vector<DAGNode *> &target,
    std::vector<DAGNode *> &internal) {
  std::map<const DAGNode *, int> retval{};
  retval[nullptr] = 0;    // we use nullptr for the count

  add_dagnode_to_index_entries(source, retval);
  add_dagnode_to_index_entries(internal, retval);
  add_dagnode_to_index_entries(target, retval);

  return retval;
}


std::vector<Edge> create_edges(std::map<const DAGNode *, int> &dtoi,
                               std::vector<FullDAGEdge> &edges) {
  std::vector<Edge> retval{};

  for (size_t i = 0; i < edges.size(); ++i) {
    retval.emplace_back(
      Edge{dtoi[edges[i].source], dtoi[edges[i].target], edges[i].op}
    );
  }

  return retval;
}


std::vector<Node> create_nodes(std::map<const DAGNode *, int> &dtoi,
                               std::vector<DAGNode *> &source,
                               std::vector<DAGNode *> &target,
                               std::vector<DAGNode *> &internal) {
  std::vector<Node> retval(dtoi.size() - 1);
  for (size_t i = 0; i < source.size(); ++i) {
    retval[dtoi[source[i]]].type = NodeType::Source;
    retval[dtoi[source[i]]].locality = source[i]->locality;
  }
  for (size_t i = 0; i < internal.size(); ++i) {
    retval[dtoi[internal[i]]].locality = internal[i]->locality;
  }
  for (size_t i = 0; i < target.size(); ++i) {
    retval[dtoi[target[i]]].type = NodeType::Target;
    retval[dtoi[target[i]]].locality = target[i]->locality;
  }
  return retval;
}


void compute_internal_types(std::map<const DAGNode *, int> &dtoi,
                            std::vector<DAGNode *> &internals,
                            std::vector<Node> &nodes) {
  // loop over the internals and see if there are edges that make it
  // clear which is which - if there are not out edges, this is likely
  // L from a BH run - actually, the source root gets incorrectly identified
  // this way
  //TODO fix that source root problem
  for (size_t i = 0; i < internals.size(); ++i) {
    if (internals[i]->edges.size()) {
      nodes[dtoi[internals[i]]].type
          = edge_type_to_node_type(internals[i]->edges[0].op);
    } else {
      nodes[dtoi[internals[i]]].type = NodeType::Local;
    }
  }
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




void output_dag_as_JSON(std::vector<DAGNode *> &source,
                        std::vector<DAGNode *> &target,
                        std::vector<DAGNode *> &internal) {
  // Create a mapping from DAGNode * to index
  auto dagnode_to_index = create_dagnode_to_index(source, target, internal);

  // Create the edge records
  auto full_edges = collect_full_dag_edges(source, target, internal);
  std::vector<Edge> E = create_edges(dagnode_to_index, full_edges);

  // Get the nodes made
  std::vector<Node> N = create_nodes(dagnode_to_index, source, target,
                                     internal);
  compute_internal_types(dagnode_to_index, internal, N);
  int n_s, n_t;
  compute_depths(N, n_s, n_t);
  create_unique_node_ids(N);

  // output
  print_json(N, E, n_s, n_t);
}


} // namespace dashmm
