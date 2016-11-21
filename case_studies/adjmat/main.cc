#include <cassert>
#include <cstdio>

#include <map>
#include <vector>


struct Index {
  int x;
  int y;
  int z;
  int level;
  int source;
};


enum class Operation {
  StoT = 0,
  StoM = 1,
  StoL = 2,
  MtoM = 3,
  MtoL = 4,
  MtoT = 5,
  MtoI = 6,
  ItoI = 7,
  ItoL = 8,
  LtoL = 9,
  LtoT = 10
};


struct Edge {
  Index source;
  Index target;
  Operation op;
  int sloc;
  int tloc;
};


struct IdTriple {
  int global;
  int source;
  int target;
};


bool operator<(const Index &a, const Index &b) {
  if (a.x < b.x) return true;
  if (a.y < b.y) return true;
  if (a.z < b.z) return true;
  if (a.level < b.level) return true;
  if (a.source < b.source) return true;
  return false;
}


// decide what is the operation
Operation string_to_op(char *code) {
  if (code[0] == 'S') {
    if (code[3] == 'T') {
      return Operation::StoT;
    } else if (code[3] == 'M') {
      return Operation::StoM;
    } else if (code[3] == 'L') {
      return Operation::StoL;
    }
  } else if (code[0] == 'M') {
    if (code[3] == 'M') {
      return Operation::MtoM;
    } else if (code[3] == 'L') {
      return Operation::MtoL;
    } else if (code[3] == 'T') {
      return Operation::MtoT;
    } else if (code[3] == 'I') {
      return Operation::MtoI;
    }
  } else if (code[0] == 'I') {
    if (code[3] == 'I') {
      return Operation::ItoI;
    } else if (code[3] == 'L') {
      return Operation::ItoL;
    }
  } else if (code[0] == 'L') {
    if (code[3] == 'L') {
      return Operation::LtoL;
    } else if (code[3] == 'T') {
      return Operation::LtoT;
    }
  }
  assert(0 && "Bad code");
  return Operation::StoT;
}


// Decide if the operation is from a source node
int source_is_source(Operation op) {
  int retval{0};
  switch (op) {
    case Operation::StoT:
    case Operation::StoM:
    case Operation::StoL:
      retval = 0;
      break;
    case Operation::MtoM:
    case Operation::MtoL:
    case Operation::MtoT:
    case Operation::MtoI:
      retval = 1;
      break;
    case Operation::ItoI:
      retval = 2;
      break;
    case Operation::ItoL:
      retval = 5;
      break;
    case Operation::LtoL:
    case Operation::LtoT:
      retval = 4;
      break;
    default:
      assert(0 && "Say what?");
      break;
  }
  return retval;
}


// Decide if the operation is to a source node
int target_is_source(Operation op) {
  int retval{0};
  switch (op) {
    case Operation::StoT:
      retval = 3;
      break;
    case Operation::StoM:
      retval = 1;
      break;
    case Operation::StoL:
      retval = 4;
      break;
    case Operation::MtoM:
      retval = 1;
      break;
    case Operation::MtoL:
      retval = 4;
      break;
    case Operation::MtoT:
      retval = 3;
      break;
    case Operation::MtoI:
      retval = 2;
      break;
    case Operation::ItoI:
      retval = 5;
      break;
    case Operation::ItoL:
      retval = 4;
      break;
    case Operation::LtoL:
      retval = 4;
      break;
    case Operation::LtoT:
      retval = 3;
      break;
    default:
      assert(0 && "No way...");
      break;
  }
  return retval;
}


// Read the edges from the file
std::vector<Edge> read_edges_from_file(char *fname) {
  std::vector<Edge> retval{};

  FILE *ifd = fopen(fname, "r");
  if (ifd == nullptr) {
    return retval;
  }

  int a, b, c, d, e, f, g, h, i, j;
  char code[20];
  while (11 == fscanf(ifd, "%d %d %d %d - %d - %s - %d %d %d %d - %d",
                      &a, &b, &c, &d, &e, code, &f, &g, &h, &i, &j)) {
    // convert string to code
    Operation op = string_to_op(code);
    retval.push_back(Edge{Index{a, b, c, d, source_is_source(op)},
                          Index{f, g, h, i, target_is_source(op)},
                          op, e, j});
  }

  fclose(ifd);

  return retval;
}


std::map<Index, IdTriple> map_index_to_id(const std::vector<Edge> &edges,
                                          int &nnodes, int &nsource,
                                          int &ntarget) {
  std::map<Index, IdTriple> retval{};
  int globalid{0};    // running id for the global list
  int sourceid{0};    // running id for the source list
  int targetid{0};    // running id for the target list

  for (auto i = edges.begin(); i != edges.end(); ++i) {
    auto selem = retval.find(i->source);
    if (selem == retval.end()) {
      auto value = IdTriple{globalid++, sourceid++, -1};
      retval[i->source] = value;
    } else if (selem->second.source == -1) {
      selem->second.source = sourceid++;
    }

    auto telem = retval.find(i->target);
    if (telem == retval.end()) {
      auto value = IdTriple{globalid++, -1, targetid++};
      retval[i->target] = value;
    } else if (telem->second.target == -1) {
      telem->second.target = targetid++;
    }
  }

  nnodes = globalid;
  nsource = sourceid;
  ntarget = targetid;

  return retval;
}


size_t xy_to_index(int s, int t, int n_target) {
  return s * n_target + t;
}


std::vector<int> create_unsorted_matrix(const std::vector<Edge> &edges,
                                        const std::map<Index, IdTriple> &ids,
                                        int n_source, int n_target) {
  std::vector<int> retval(n_source * n_target, 0);

  for (auto i = edges.begin(); i != edges.end(); ++i) {
    // find source and target ids
    auto selem = ids.find(i->source);
    assert(selem != ids.end());
    auto telem = ids.find(i->target);
    assert(telem != ids.end());

    // set appropriate element of retval
    size_t idx = xy_to_index(selem->second.source, telem->second.target,
                             n_target);
    retval[idx] = 1;
  }

  return retval;
}


void print_unsorted_matrix(const std::vector<int> &mat,
                           int n_source, int n_target, char *fname) {
  FILE *ofd = fopen(fname, "w");
  assert(ofd != nullptr);

  size_t offset = 0;
  for (int i = 0; i < n_source; ++i) {
    for (int j = 0; j < n_target; ++j) {
      fprintf(ofd, "%d ", mat[offset]);
      offset++;
    }
  }

  fclose(ofd);
}


int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stdout, "usage: %s <input csv edge file> <out csv>\n", argv[0]);
    return 0;
  }

  // load the file
  auto edges = read_edges_from_file(argv[1]);

  // map indices into unique ids
  int n_total{0};
  int n_source{0};
  int n_target{0};
  auto idxmap = map_index_to_id(edges, n_total, n_source, n_target);

  // print out the matrix
  auto unsort = create_unsorted_matrix(edges, idxmap, n_source, n_target);
  print_unsorted_matrix(unsort, n_source, n_target, argv[2]);

  return 0;
}