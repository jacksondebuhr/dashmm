#include <cassert>
#include <cstdio>

#include <algorithm>
#include <map>
#include <vector>


struct Triple {
  int gid;
  int sid;
  int tid;

  int loc;
  int sorted;

  Triple() : gid{-1}, sid{-1}, tid{-1}, loc{-1}, sorted{-1} { }
};

struct Edge {
  uint64_t source;
  uint64_t target;
  int sloc;
  int tloc;
  int op;
};

struct NodeLocality {
  uint64_t node;
  int loc;
};


bool compare_node_locality(const NodeLocality &a, const NodeLocality &b) {
  return a.loc < b.loc;
}


int op_from_code(char *code) {
  if (code[0] == 'S') {
    if (code[3] == 'T') {
      return 1;
    } else if (code[3] == 'M') {
      return 2;
    } else if (code[3] == 'L') {
      return 3;
    }
  } else if (code[0] == 'M') {
    if (code[3] == 'M') {
      return 4;
    } else if (code[3] == 'L') {
      return 5;
    } else if (code[3] == 'T') {
      return 6;
    } else if (code[3] == 'I') {
      return 7;
    }
  } else if (code[0] == 'I') {
    if (code[3] == 'I') {
      return 8;
    } else if (code[3] == 'L') {
      return 9;
    }
  } else if (code[0] == 'L') {
    if (code[3] == 'L') {
      return 10;
    } else if (code[3] == 'T') {
      return 11;
    }
  }
  assert(0 && "Bad code");
  return 0;
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
  uint64_t sptr, tptr;
  while (13 == fscanf(ifd, "%d %d %d %d - %d - %s - %d %d %d %d - %d - %lx %lx",
                      &a, &b, &c, &d, &e, code, &f, &g, &h, &i, &j,
                      &sptr, &tptr)) {
    // convert string to code
    retval.push_back(Edge{sptr, tptr, e, j, op_from_code(code)});
  }

  fclose(ifd);

  return retval;
}

// map from node address to indices
std::map<uint64_t, Triple> map_node_to_id(const std::vector<Edge> &edges,
                                          int &n_total, int &n_source,
                                          int &n_target) {
  std::map<uint64_t, Triple> retval{};

  // assign global index
  for (size_t i = 0; i < edges.size(); ++i) {
    auto selem = retval.find(edges[i].source);
    if (selem == retval.end()) {
      retval[edges[i].source].gid = n_total++;
      retval[edges[i].source].loc = edges[i].sloc;
    }

    auto telem = retval.find(edges[i].target);
    if (telem == retval.end()) {
      retval[edges[i].target].gid = n_total++;
      retval[edges[i].target].loc = edges[i].tloc;
    }
  }

  // assign source and target
  for (size_t i = 0; i < edges.size(); ++i) {
    assert(retval.count(edges[i].source));
    if (retval[edges[i].source].sid == -1) {
      retval[edges[i].source].sid = n_source++;
    }
    assert(retval.count(edges[i].target));
    if (retval[edges[i].target].tid == -1) {
      retval[edges[i].target].tid = n_target++;
    }
  }

  return retval;
}

size_t xy_to_index(int s, int t, int n_target) {
  return s * n_target + t;
}


std::vector<int> create_sourcetarget_matrix(const std::vector<Edge> &edges,
                                        const std::map<uint64_t, Triple> &ids,
                                        int n_source, int n_target) {
  std::vector<int> retval(n_source * n_target, 0);

  for (size_t i = 0; i < edges.size(); ++i) {
    // find source and target ids
    auto selem = ids.find(edges[i].source);
    assert(!(selem == ids.end()));
    auto sidx = selem->second.sid;

    auto telem = ids.find(edges[i].target);
    assert(!(telem == ids.end()));
    auto tidx = telem->second.tid;

    // set appropriate element of retval
    size_t idx = xy_to_index(sidx, tidx, n_target);
    retval[idx] = edges[i].op;
  }

  return retval;
}


std::vector<int> create_unsorted_matrix(const std::vector<Edge> &edges,
                                        const std::map<uint64_t, Triple> &ids,
                                        int n_total) {
  std::vector<int> retval(n_total * n_total, 0);

  for (size_t i = 0; i < edges.size(); ++i) {
    // find source and target ids
    auto selem = ids.find(edges[i].source);
    assert(!(selem == ids.end()));
    auto sidx = selem->second.gid;

    auto telem = ids.find(edges[i].target);
    assert(!(telem == ids.end()));
    auto tidx = telem->second.gid;

    // set appropriate element of retval
    size_t idx = xy_to_index(sidx, tidx, n_total);
    retval[idx] = edges[i].op;
  }

  return retval;
}

std::vector<int> create_sorted_matrix(const std::vector<Edge> &edges,
                                      const std::map<uint64_t, Triple> &ids,
                                      int n_total) {
  std::vector<int> retval(n_total * n_total, 0);

  for (size_t i = 0; i < edges.size(); ++i) {
    // find source and target ids
    auto selem = ids.find(edges[i].source);
    assert(!(selem == ids.end()));
    auto sidx = selem->second.sorted;

    auto telem = ids.find(edges[i].target);
    assert(!(telem == ids.end()));
    auto tidx = telem->second.sorted;

    // set appropriate element of retval
    size_t idx = xy_to_index(sidx, tidx, n_total);
    assert(edges[i].tloc > -1);
    assert(edges[i].sloc > -1);
    retval[idx] = edges[i].tloc + 6 * edges[i].sloc + 1;
  }

  return retval;
}


void print_matrix(const std::vector<int> &mat,
                  int n_source, int n_target, char *fname, char *extend) {
  char buffer[300];
  sprintf(buffer, "%s-%s", fname, extend);

  FILE *ofd = fopen(buffer, "w");
  assert(ofd != nullptr);

  size_t offset = 0;
  for (int i = 0; i < n_source; ++i) {
    for (int j = 0; j < n_target; ++j) {
      fprintf(ofd, "%d ", mat[offset]);
      offset++;
    }
    fprintf(ofd, "\n");
  }

  fclose(ofd);
}


void fill_in_sorted_ids(std::map<uint64_t, Triple> &ids) {
  std::vector<NodeLocality> shuffle{};

  // put pairs into a vector
  for (auto i = ids.begin(); i != ids.end(); ++i) {
    shuffle.emplace_back(NodeLocality{i->first, i->second.loc});
  }

  // sort by locality
  std::sort(shuffle.begin(), shuffle.end(), compare_node_locality);

  // assign sorted ids to the original map
  int sortedid{0};
  for (size_t i = 0; i < shuffle.size(); ++i) {
    auto elem = ids.find(shuffle[i].node);
    assert(elem != ids.end());
    elem->second.sorted = sortedid++;
  }
}


int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stdout, "usage: %s <input csv edge file> <out csv base>\n", argv[0]);
    return 0;
  }

  // load the file
  auto edges = read_edges_from_file(argv[1]);

  // map indices into unique ids
  int n_total{0};
  int n_source{0};
  int n_target{0};
  auto idx = map_node_to_id(edges, n_total, n_source, n_target);

  // print out the matrix arranged by source and target
  auto stmat = create_sourcetarget_matrix(edges, idx, n_source, n_target);
  print_matrix(stmat, n_source, n_target, argv[2], "st.txt");

  // print out the unsorted global matrix
  auto unsort = create_unsorted_matrix(edges, idx, n_total);
  print_matrix(unsort, n_total, n_total, argv[2], "g.txt");

  // Create the node and locality map
  fill_in_sorted_ids(idx);

  // print out the sorted global matrix
  auto sortmat = create_sorted_matrix(edges, idx, n_total);
  print_matrix(sortmat, n_total, n_total, argv[2], "sort.txt");

  return 0;
}