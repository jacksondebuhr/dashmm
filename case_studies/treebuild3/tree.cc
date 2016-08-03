#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <functional>
#include "tree.h"

double corner_x;        /// The corner of the overall domain in x
double corner_y;        /// -in y
double corner_z;        /// -in z
double size;            /// The size of the overall domain

// TODO: Many of these are implemented as cyclic arrays that point to the local
// data. But is there any reason to do this? We never make use of the fact that
// these are in the global address space. All the action is local. There may be
// exceptions on this, but it is worth looking at. Unless we actually go ahead
// and use the GAS nature of this, we are just making some things harder on
// ourselves.
int unif_level;         /// The level of uniform partition
hpx_addr_t unif_count;  /// A cyclic array for counting uniform parition
                        /// occupation
hpx_addr_t unif_grid;   /// The uniform grid; that is, a cyclic array of
                        /// MetaData pointing to local arrays of Nodes.
hpx_addr_t unif_done;   /// A cyclic array for detecting completion of the
                        /// uniform portion of tree creation

hpx_addr_t sorted_src;  /// Array meta data for the sorted sources
hpx_addr_t sorted_tar;  /// Array meta data for the sorted targets
int *distribute;        /// Array giving the distribution of the uniform grid

int threshold;
int *swap_src;
int *bin_src;
int *map_src;
int *swap_tar;
int *bin_tar;
int *map_tar;

void set_point_in_cube(Point &p) {
  double x = 1.0 * rand() / RAND_MAX - 0.5;
  double y = 1.0 * rand() / RAND_MAX - 0.5;
  double z = 1.0 * rand() / RAND_MAX - 0.5;
  p = Point{x, y, z};
}

void set_point_on_sphere(Point &p) {
  double theta = 1.0 * rand() / RAND_MAX * M_PI_2;
  double phi = 1.0 * rand() / RAND_MAX * M_PI * 2;
  double x = sin(theta) * cos(phi);
  double y = sin(theta) * sin(phi);
  double z = cos(theta);
  p = Point{x, y, z};
}

uint64_t split(unsigned k) {
  uint64_t split = k & 0x1fffff;
  split = (split | split << 32) & 0x1f00000000ffff;
  split = (split | split << 16) & 0x1f0000ff0000ff;
  split = (split | split << 8)  & 0x100f00f00f00f00f;
  split = (split | split << 4)  & 0x10c30c30c30c30c3;
  split = (split | split << 2)  & 0x1249249249249249;
  return split;
}

uint64_t morton_key(unsigned x, unsigned y, unsigned z) {
  uint64_t key = 0;
  key |= split(x) | split(y) << 1 | split(z) << 2;
  return key;
}

// Given the global counts, this will partition the uniform grid among the
// available localities. There is perhaps some room for simplification here.
// I would rather shoot for fixed targets instead of aiming to take a fair
// fraction of whatever is left. But there might not be any real performance
// impact either way.
//
// TODO: Note that this would be a target for a second member of the
// Distribution Policy
int *distribute_points(int num_ranks, const int *global, int len) {
  int *ret = new int[num_ranks]();

  const int *s = global; // Source counts
  const int *t = &global[len]; // Target counts

  int total = 0;
  for (int i = 0; i < len; ++i) {
    total += s[i] + t[i];
  }

  int rank = 0;
  int iterator = 0;

  while (rank < num_ranks && total > 0) {
    if (rank == num_ranks - 1) {
      // Take the remaining grids
      ret[rank++] = len - 1;
    } else {
      int avg = total / (num_ranks - rank);
      int sum = 0;
      int sum1 = 0;

      for (int i = iterator; i < len; ++i) {
        sum += s[i] + t[i];
        if (i == len - 1) {
          // There will be ranks left without grids assigned.
          delete [] ret;
          return nullptr;
        } else {
          sum1 = sum + s[i + 1] + t[i + 1];
        }

        if (sum <= avg && avg <= sum1) {
          // Check which is closer to avg
          if (avg - sum <= sum1 - avg) {
            iterator = i;
          } else {
            iterator = i + 1;
            sum = sum1;
          }
          break;
        }
      }

      ret[rank++] = iterator;
      total -= sum;
      iterator++;
    }
  }

  return ret;
}

void generate_weak_scaling_input(ArrayMetaData *meta_s, int nsrc_per_rank,
                                 ArrayMetaData *meta_t, int ntar_per_rank,
                                 char datatype, int seed) {
  meta_s->size = sizeof(Point);
  meta_s->count = nsrc_per_rank;
  meta_s->data = new char[sizeof(Point) * nsrc_per_rank]();

  meta_t->size = sizeof(Point);
  meta_t->count = ntar_per_rank;
  meta_t->data = new char[sizeof(Point) * ntar_per_rank]();

  Point *p_s = reinterpret_cast<Point *>(meta_s->data);
  Point *p_t = reinterpret_cast<Point *>(meta_t->data);

  srand(seed);
  std::function<void(Point &)> set_point =
    (datatype == 'c' ? set_point_in_cube : set_point_on_sphere);

  for (int i = 0; i < nsrc_per_rank; ++i)
    set_point(p_s[i]);

  for (int i = 0; i < ntar_per_rank; ++i)
    set_point(p_t[i]);
}

void generate_strong_scaling_input(ArrayMetaData *meta_s, int nsrc,
                                   ArrayMetaData *meta_t, int ntar,
                                   char datatype, int rank, int num_ranks,
                                   int nseed) {
  // Divide nsrc and ntar into nseed chunks. Use seed (2 + r) to populate
  // the points in each chunk, where 0 <= r < nseed

  std::function<void(Point &)> set_point =
    (datatype == 'c' ? set_point_in_cube : set_point_on_sphere);

  int q1, r1, q2, r2, nsrc_curr_rank, s1, s2, ntar_curr_rank, t1, t2;

  // Generate source points
  // nsrc = q1 * num_ranks + r1 = q2 * nseed + r2
  q1 = nsrc / num_ranks;
  r1 = nsrc % num_ranks;
  q2 = nsrc / nseed;
  r2 = nsrc % nseed;
  nsrc_curr_rank = (rank < r1 ? q1 + 1 : q1);
  s1 = (rank < r1 ? (q1 + 1) * rank : (q1 + 1) * r1 + q1 * (rank - r1));
  s2 = s1 + nsrc_curr_rank - 1;

  meta_s->size = sizeof(Point);
  meta_s->count = nsrc_curr_rank;
  meta_s->data = new char[sizeof(Point) * nsrc_curr_rank]();
  Point *p_s = reinterpret_cast<Point *>(meta_s->data);

  for (int r = 0; r < nseed; ++r) {
    int nsrc_curr_seed = (r < r2 ? q2 + 1 : q2);
    int s3 = (r < r2 ? (q2 + 1) * r : (q2 + 1) * r2 + q2 * (r - r2));
    int s4 = s3 + nsrc_curr_seed - 1;
    if (s4 < s1) {
      continue;
    } else if (s3 > s2) {
      break;
    } else {
      srand(r + 2);
      for (int s = s3; s <= s4; ++s) {
        if (s >= s1 && s <= s2)
          set_point(p_s[s - s1]);
      }
    }
  }

  // Generate target points
  // ntar = q1 * num_ranks + r1 = q2 * nseed + r2
  q1 = ntar / num_ranks;
  r1 = ntar % num_ranks;
  q2 = ntar / nseed;
  r2 = ntar % nseed;
  ntar_curr_rank = (rank < r1 ? q1 + 1 : q1);
  t1 = (rank < r1 ? (q1 + 1) * rank : (q1 + 1) * r1 + q1 * (rank - r1));
  t2 = t1 + ntar_curr_rank - 1;

  meta_t->size = sizeof(Point);
  meta_t->count = ntar_curr_rank;
  meta_t->data = new char[sizeof(Point) * ntar_curr_rank]();
  Point *p_t = reinterpret_cast<Point *>(meta_t->data);

  for (int r = 0; r < nseed; ++r) {
    int ntar_curr_seed = (r < r2 ? q2 + 1 : q2);
    int t3 = (r < r2 ? (q2 + 1) * r : (q2 + 1) * r2 + q2 * (r - r2));
    int t4 = t3 + ntar_curr_seed - 1;
    if (t4 < t1) {
      continue;
    } else if (t3 > t2) {
      break;
    } else {
      srand(r + 2);
      for (int t = t3; t <= t4; ++t) {
        if (t >= t1 && t <= t2)
          set_point(p_t[t - t1]);
      }
    }
  }
}

/// This will generate the points for the trees.
/// This is pretty simple. The local offset to the correct part of the
/// meta data for the array is computed, and the values are filled in.
/// nothing fancy.
int generate_input_handler(char scaling, char datatype, int nsrc, int ntar,
                           int nseed, hpx_addr_t sources, hpx_addr_t targets) {
  int rank = hpx_get_my_rank();
  int num_ranks = hpx_get_num_ranks();
  hpx_addr_t curr_s = hpx_addr_add(sources, sizeof(ArrayMetaData) * rank,
                                   sizeof(ArrayMetaData));
  hpx_addr_t curr_t = hpx_addr_add(targets, sizeof(ArrayMetaData) * rank,
                                   sizeof(ArrayMetaData));

  ArrayMetaData *meta_s{nullptr}, *meta_t{nullptr};
  if (!hpx_gas_try_pin(curr_s, (void **)&meta_s) ||
      !hpx_gas_try_pin(curr_t, (void **)&meta_t))
    return HPX_ERROR;

  if (scaling == 'w') {
    // Generate input for weak scaling test
    // nsrc and ntar are the number of sources and targets per rank

    // Set the seed of the random number generator. Note that if the seed is set
    // to 1, the generator is reinitialized to its initial value. For this
    // reason, the seed is chosen to be 2 + rank here.
    generate_weak_scaling_input(meta_s, nsrc, meta_t, ntar, datatype, rank + 2);
  } else {
    // Generate input for strong scaling test
    // nsrc and ntar are the total number of sources and targets
    generate_strong_scaling_input(meta_s, nsrc, meta_t, ntar, datatype,
                                  rank, num_ranks, nseed);
  }
  hpx_gas_unpin(curr_s);
  hpx_gas_unpin(curr_t);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, generate_input_action, generate_input_handler,
           HPX_CHAR, HPX_CHAR, HPX_INT, HPX_INT, HPX_INT, HPX_ADDR, HPX_ADDR);



/// Pretty simple really, just gets rid of the sources and targets
int destroy_input_handler(hpx_addr_t sources, hpx_addr_t targets) {
  int rank = hpx_get_my_rank();
  hpx_addr_t curr_s = hpx_addr_add(sources, sizeof(ArrayMetaData) * rank,
                                   sizeof(ArrayMetaData));
  hpx_addr_t curr_t = hpx_addr_add(targets, sizeof(ArrayMetaData) * rank,
                                   sizeof(ArrayMetaData));

  ArrayMetaData *meta_s{nullptr}, *meta_t{nullptr};
  if (!hpx_gas_try_pin(curr_s, (void **)&meta_s) ||
      !hpx_gas_try_pin(curr_t, (void **)&meta_t))
    return HPX_ERROR;

  delete [] meta_s->data;
  delete [] meta_t->data;

  hpx_gas_unpin(curr_s);
  hpx_gas_unpin(curr_t);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, destroy_input_action, destroy_input_handler,
           HPX_ADDR, HPX_ADDR);

void domain_geometry_init_handler(void *data, const size_t size,
                                  int *init, size_t init_size) {
  char *input = static_cast<char *>(data);
  int *cnt = reinterpret_cast<int *>(input);
  double *values = reinterpret_cast<double *>(input + sizeof(int));

  *cnt = *init; // # of inputs

  values[0] = 1e50; // xmin
  values[1] = -1e50; // xmax
  values[2] = 1e50; // ymin
  values[3] = -1e50; // ymax
  values[4] = 1e50; // zmin
  values[5] = -1e50; // zmax
}
HPX_ACTION(HPX_FUNCTION, 0, domain_geometry_init_action,
           domain_geometry_init_handler, HPX_POINTER, HPX_SIZE_T,
           HPX_POINTER, HPX_SIZE_T);

void domain_geometry_op_handler(void *data, double *rhs, size_t size) {
  char *input = static_cast<char *>(data);
  int *cnt = reinterpret_cast<int *>(input);
  double *lhs = reinterpret_cast<double *>(input + sizeof(int));

  (*cnt)--;
  lhs[0] = fmin(lhs[0], rhs[0]);
  lhs[1] = fmax(lhs[1], rhs[1]);
  lhs[2] = fmin(lhs[2], rhs[2]);
  lhs[3] = fmax(lhs[3], rhs[3]);
  lhs[4] = fmin(lhs[4], rhs[4]);
  lhs[5] = fmax(lhs[5], rhs[5]);
}
HPX_ACTION(HPX_FUNCTION, 0, domain_geometry_op_action,
           domain_geometry_op_handler, HPX_POINTER, HPX_POINTER, HPX_SIZE_T);

bool domain_geometry_predicate_handler(void *data, size_t size) {
  char *input = static_cast<char *>(data);
  int *cnt = reinterpret_cast<int *>(input);
  return (*cnt == 0);
}
 HPX_ACTION(HPX_FUNCTION, 0, domain_geometry_predicate_action,
            domain_geometry_predicate_handler, HPX_POINTER, HPX_SIZE_T);

int set_domain_geometry_handler(hpx_addr_t sources,
                                hpx_addr_t targets,
                                hpx_addr_t domain_geometry) {
  int rank = hpx_get_my_rank();
  hpx_addr_t curr_s = hpx_addr_add(sources, sizeof(ArrayMetaData) * rank,
                                   sizeof(ArrayMetaData));
  hpx_addr_t curr_t = hpx_addr_add(targets, sizeof(ArrayMetaData) * rank,
                                   sizeof(ArrayMetaData));
  ArrayMetaData *meta_s{nullptr}, *meta_t{nullptr};

  if (!hpx_gas_try_pin(curr_s, (void **)&meta_s) ||
      !hpx_gas_try_pin(curr_t, (void **)&meta_t))
    return HPX_ERROR;

  Point *s = reinterpret_cast<Point *>(meta_s->data);
  Point *t = reinterpret_cast<Point *>(meta_t->data);
  int nsrc = meta_s->count;
  int ntar = meta_t->count;

  double var[6] = {1e50, -1e50, 1e50, -1e50, 1e50, -1e50};

  // NOTE: Here is an opportunity to do even more parallelism. This will be
  // a single thread per locality. Why not do the on locality reduction in
  // parallel.
  for (int i = 0; i < nsrc; ++i) {
    var[0] = fmin(var[0], s[i].x());
    var[1] = fmax(var[1], s[i].x());
    var[2] = fmin(var[2], s[i].y());
    var[3] = fmax(var[3], s[i].y());
    var[4] = fmin(var[4], s[i].z());
    var[5] = fmax(var[5], s[i].z());
  }

  for (int i = 0; i < ntar; ++i) {
    var[0] = fmin(var[0], t[i].x());
    var[1] = fmax(var[1], t[i].x());
    var[2] = fmin(var[2], t[i].y());
    var[3] = fmax(var[3], t[i].y());
    var[4] = fmin(var[4], t[i].z());
    var[5] = fmax(var[5], t[i].z());
  }

  hpx_gas_unpin(curr_s);
  hpx_gas_unpin(curr_t);
  hpx_lco_set_lsync(domain_geometry, sizeof(double) * 6, var, HPX_NULL);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, set_domain_geometry_action,
           set_domain_geometry_handler, HPX_ADDR, HPX_ADDR, HPX_ADDR);

void unif_count_init_handler(void *data, const size_t size,
                             int *init, size_t init_size) {

  int *input = reinterpret_cast<int *>(static_cast<char *>(data));

  input[0] = init[0]; // # of inputs;
  input[1] = init[1]; // # of terms

  for (int i = 0; i < init[1]; ++i)
    input[i + 2] = 0;
}
HPX_ACTION(HPX_FUNCTION, 0, unif_count_init_action,
           unif_count_init_handler, HPX_POINTER, HPX_SIZE_T,
           HPX_POINTER, HPX_SIZE_T);

void unif_count_op_handler(void *data, int *rhs, size_t size) {
  int *input = reinterpret_cast<int *>(static_cast<char *>(data));
  int *lhs = &input[2];

  input[0]--;
  int nterms = input[1];
  for (int i = 0; i < nterms; ++i)
    lhs[i] += rhs[i];
}
HPX_ACTION(HPX_FUNCTION, 0, unif_count_op_action,
           unif_count_op_handler, HPX_POINTER, HPX_POINTER, HPX_SIZE_T);

bool unif_count_predicate_handler(void *data, size_t size) {
  int *input = reinterpret_cast<int *>(static_cast<char *>(data));
  return (input[0] == 0);
}
HPX_ACTION(HPX_FUNCTION, 0, unif_count_predicate_action,
           unif_count_predicate_handler, HPX_POINTER, HPX_SIZE_T);

int init_partition_handler(hpx_addr_t count, hpx_addr_t done, hpx_addr_t grid,
                           hpx_addr_t sources, hpx_addr_t targets, int level,
                           int limit, double cx, double cy, double cz,
                           double sz) {
  int rank = hpx_get_my_rank();
  int num_ranks = hpx_get_num_ranks();
  int dim = pow(2, level), dim3 = pow(8, level);
  threshold = limit;

  // Here we save the addresses of these object at the other ranks.
  if (rank) {
    unif_count = count;
    unif_done = done;
    unif_grid = grid;
    sorted_src = sources;
    sorted_tar = targets;
    unif_level = level;
    corner_x = cx;
    corner_y = cy;
    corner_z = cz;
    size = sz;
  }

  // TODO: currently, this unif_count is an array of such LCOs. Is there
  // some reason to think that this is better than just a single LCO that
  // everyone knows the address of? Perhaps there is slightly more asynchrony
  // with this. But there is more information to send around.
  //
  // TODO: Also, there is no need for this to be a user LCO; the built in
  // reduction LCO would be just fine here.

  // Setup unif_count LCO
  hpx_addr_t curr_unif_count = hpx_addr_add(unif_count,
                                            sizeof(hpx_addr_t) * rank,
                                            sizeof(hpx_addr_t));
  hpx_addr_t *user_lco{nullptr};
  if (!hpx_gas_try_pin(curr_unif_count, (void **)&user_lco))
    return HPX_ERROR;

  int init[2] = {num_ranks, dim3 * 2};
  *user_lco = hpx_lco_user_new(sizeof(int) * (2 + dim3 * 2),
                               unif_count_init_action,
                               unif_count_op_action,
                               unif_count_predicate_action,
                               &init, sizeof(init));
  hpx_gas_unpin(curr_unif_count);

  // NOTE: Unlike the previous, this one should be on a rank-by-rank basis,
  // so that future stuff can begin asap.

  // Setup unif_done LCO
  hpx_addr_t curr_unif_done = hpx_addr_add(unif_done,
                                           sizeof(hpx_addr_t) * rank,
                                           sizeof(hpx_addr_t));
  hpx_addr_t *gate{nullptr};
  if (!hpx_gas_try_pin(curr_unif_done, (void **)&gate))
    return HPX_ERROR;

  *gate = hpx_lco_and_new(1);
  hpx_gas_unpin(curr_unif_done);

  // NOTE: The uniform grid means an array of the nodes that give the
  // uniform partition. Each rank will have the same nodes in the array.
  // NOTE: Is there a good reason for these to be part of an array in GAS?
  // Is that important in some way?

  // Setup unif_grid
  hpx_addr_t curr_unif_grid = hpx_addr_add(unif_grid,
                                           sizeof(ArrayMetaData) * rank,
                                           sizeof(ArrayMetaData));
  ArrayMetaData *meta{nullptr};
  if (!hpx_gas_try_pin(curr_unif_grid, (void **)&meta))
    return HPX_ERROR;

  meta->size = sizeof(Node);
  meta->count = 2 * dim3;
  meta->data = new char[sizeof(Node) * 2 * dim3]();

  Node *n = reinterpret_cast<Node *>(meta->data);
  for (int iz = 0; iz < dim; ++iz) {
    for (int iy = 0; iy < dim; ++iy) {
      for (int ix = 0; ix < dim; ++ix) {
        uint64_t mid = morton_key(ix, iy, iz);
        // TODO: Are these copies needed? Why not just set the index?
        n[mid] = Node{Index{level, ix, iy, iz}};
        n[mid + dim3] = Node{Index{level, ix, iy, iz}};
      }
    }
  }

  hpx_gas_unpin(curr_unif_grid);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, init_partition_action, init_partition_handler,
           HPX_ADDR, HPX_ADDR, HPX_ADDR, HPX_ADDR, HPX_ADDR, HPX_INT, HPX_INT,
           HPX_DOUBLE, HPX_DOUBLE, HPX_DOUBLE, HPX_DOUBLE);

// This is the target action of the count exchange. This is called n-ranks
// times, resulting in n-ranks - 1 network messages for each.
int exchange_count_handler(void *args, size_t size) {
  int rank = hpx_get_my_rank();
  hpx_addr_t curr_unif_count = hpx_addr_add(unif_count,
                                            sizeof(hpx_addr_t) * rank,
                                            sizeof(hpx_addr_t));
  hpx_addr_t *user_lco{nullptr};
  if (!hpx_gas_try_pin(curr_unif_count, (void **)&user_lco)) {
    return HPX_ERROR;
  }

  hpx_lco_set(*user_lco, size, args, HPX_NULL, HPX_NULL);
  hpx_gas_unpin(curr_unif_count);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED, exchange_count_action,
           exchange_count_handler, HPX_POINTER, HPX_SIZE_T);

// This will assign the points to the uniform grid. This gives the points
// the id (in the Morton Key sense) of the box to which they are assigned,
// and it will count the numbers in each box.
void assign_points_to_unif_grid(const Point *P, int npts, int *gid,
                                int *count, double scale) {
  // TODO: This is serial processing; is there some way to parallelize this?
  //   This would perhaps be worth timing.
  int dim = pow(2, unif_level);
  for (int i = 0; i < npts; ++i) {
    const Point *p = &P[i];
    int xid = std::min(dim - 1, (int)(dim * (p->x() - corner_x) * scale));
    int yid = std::min(dim - 1, (int)(dim * (p->y() - corner_y) * scale));
    int zid = std::min(dim - 1, (int)(dim * (p->z() - corner_z) * scale));
    gid[i] = morton_key(xid, yid, zid);
    count[gid[i]]++;
  }
}

// This will rearrange the particles into their bin order. This is a stable
// reordering.
//
// TODO: This is a sort of operation that occurs multiple times here. We should
// abstract the notion of the bin-sort (or perhaps this is bucket-sort, but
// whatever) so that we can improve the performance of that more easily as
// needed.
//
// TODO: Note that this too is a serial operation. Could we do some on-rank
// parallelism here?
int *group_points_on_unif_grid(const Point *p_in, int npts,
                               const int *gid_of_points, const int *count,
                               Point *p_out) {
  int dim3 = pow(8, unif_level);
  int *offset = new int[dim3]();
  int *assigned = new int[dim3]();

  for (int i = 1; i < dim3; ++i) {
    offset[i] = offset[i - 1] + count[i - 1];
  }

  // TODO: I think we can do away with assigned here. We can just update the
  // offset for the given bin.
  for (int i = 0; i < npts; ++i) {
    int gid = gid_of_points[i];
    int cursor = offset[gid] + assigned[gid];
    p_out[cursor] = p_in[i];
    assigned[gid]++;
  }

  delete [] assigned;
  return offset;
}

/// This is a thin wrapper around the partition method of a node. This is
/// made into an action to make use of the parallelism available with HPX.
int partition_node_handler(Node *n, Point *p, int *swap, int *bin, int *map) {
  n->partition(p, swap, bin, map, threshold, corner_x, corner_y,
               corner_z, size);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, partition_node_action, partition_node_handler,
           HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_POINTER);

// This sets up some orgnaizational structures. The return value is an
// array of offsets in the final set of points for each part of the uniform
// grid. This is basically just setup. There is a chance that some work occurs
// in one branch. Before leaving, this will bring the points that do not have
// to change rank into their correct location. If it is detected that all of
// the points for that part of the tree have arrived, then the partitioning
// work will begin.
int *init_point_exchange(int rank, const int *global_count,
                         const int *local_count, const int *local_offset,
                         const Point *temp, Node *n, ArrayMetaData *meta,
                         char type) {
  int first = (rank == 0 ? 0 : distribute[rank - 1] + 1);
  int last = distribute[rank];
  int range = last - first + 1;
  int **swap = (type == 's' ? &swap_src : &swap_tar);
  int **map = (type == 's' ? &map_src : &map_tar);
  int **bin = (type == 's' ? &bin_src : &bin_tar);

  // Compute global_offset
  int *global_offset = new int[range]();
  int num_points = global_count[first];
  for (int i = first + 1; i <= last; ++i) {
    num_points += global_count[i];
    global_offset[i - first] = global_offset[i - first - 1] +
      global_count[i - 1];
  }

  meta->size = sizeof(Point);
  meta->count = num_points;

  if (num_points > 0) {
    meta->data = new char[sizeof(Point) * num_points]();
    Point *p = reinterpret_cast<Point *>(meta->data);

    *swap = new int[num_points]();
    *map = new int[num_points]();
    *bin = new int[num_points]();
    for (int i = 0; i < num_points; ++i) {
      (*map)[i] = i;
    }

    for (int i = first; i <= last; ++i) {
      Node *curr = &n[i];
      if (global_count[i]) {
        curr->set_first(global_offset[i - first]);
        if (i < last) {
          curr->set_last(global_offset[i + 1 - first] - 1);
        } else {
          curr->set_last(num_points - 1);
        }
      }

      if (local_count[i]) {
        // Copy local points before merging remote points
        memcpy(p + global_offset[i - first], &temp[local_offset[i]],
               sizeof(Point) * local_count[i]);

        if (local_count[i] == global_count[i]) {
          // This grid does not expect remote points.
          // Spawn adaptive partitioning
          hpx_call(HPX_HERE, partition_node_action, HPX_NULL,
                   &curr, &p, swap, bin, map);
        } else {
          // Now use first_ as an iterator to trace the position to merge
          // remote points
          curr->set_first(global_offset[i - first] + local_count[i]);
        }
      }
    }
  } else {
    meta->data = nullptr;
  }

  return global_offset;
}

// This action merges particular points with the sorted list. Also, if this
// is the last set of points that are merged, this will go ahead and start
// the adaptive partitioning of that part of the local tree.
int merge_points_handler(Point *p, Point *temp, Node *n,
                         int n_arrived, int n_total, char type) {
  // Note: all the pointers are local to the calling rank.
  hpx_lco_sema_p(n->sema());
  int first = n->first();
  memcpy(p + first, temp, sizeof(Point) * n_arrived);
  n->set_first(first + n_arrived);

  if (first + n_arrived > n->last()) {
    // Spawn adaptive partitioning
    n->set_first(n->last() + 1 - n_total);

    if (type == 's') {
      hpx_call(HPX_HERE, partition_node_action, HPX_NULL,
               &n, &p, &swap_src, &bin_src, &map_src);
    } else {
      hpx_call(HPX_HERE, partition_node_action, HPX_NULL,
               &n, &p, &swap_tar, &bin_tar, &map_tar);
    }
  }
  hpx_lco_sema_v(n->sema(), HPX_NULL);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, merge_points_action, merge_points_handler,
           HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_INT, HPX_INT, HPX_CHAR);

// This is the 'far-side' of the send points message. This action merges the
// incoming points into the sorted list and will then spawn the adaptive
// partition if this happens to be the last block for a given uniform grid.
int recv_points_handler(void *args, size_t size) {
  int rank = hpx_get_my_rank();
  hpx_addr_t curr_unif_done = hpx_addr_add(unif_done,
                                           sizeof(hpx_addr_t) * rank,
                                           sizeof(hpx_addr_t));
  hpx_addr_t *gate{nullptr};
  if (!hpx_gas_try_pin(curr_unif_done, (void **)&gate))
    return HPX_ERROR;

  // Wait until the buffer is allocated before merging incoming messages
  // TODO: Is there perhaps a way to do this action as a call-when on this gate?
  hpx_lco_wait(*gate);
  hpx_gas_unpin(curr_unif_done);

  // Now merge incoming message
  hpx_addr_t curr_unif_count = hpx_addr_add(unif_count,
                                            sizeof(hpx_addr_t) * rank,
                                            sizeof(hpx_addr_t));
  hpx_addr_t curr_sorted_src = hpx_addr_add(sorted_src,
                                            sizeof(ArrayMetaData) * rank,
                                            sizeof(ArrayMetaData));
  hpx_addr_t curr_sorted_tar = hpx_addr_add(sorted_tar,
                                            sizeof(ArrayMetaData) * rank,
                                            sizeof(ArrayMetaData));
  hpx_addr_t curr_unif_grid = hpx_addr_add(unif_grid,
                                           sizeof(ArrayMetaData) * rank,
                                           sizeof(ArrayMetaData));
  hpx_addr_t *user_lco{nullptr};
  ArrayMetaData *meta_s{nullptr}, *meta_t{nullptr}, *meta_g{nullptr};
  if (!hpx_gas_try_pin(curr_unif_count, (void **)&user_lco) ||
      !hpx_gas_try_pin(curr_sorted_src, (void **)&meta_s) ||
      !hpx_gas_try_pin(curr_sorted_tar, (void **)&meta_t) ||
      !hpx_gas_try_pin(curr_unif_grid, (void **)&meta_g))
    return HPX_ERROR;

  void *user_lco_buffer{nullptr};
  hpx_lco_getref(*user_lco, 1, &user_lco_buffer);
  int *count = reinterpret_cast<int *>(static_cast<char *>(user_lco_buffer)
                                       + sizeof(int) * 2);

  Point *s = reinterpret_cast<Point *>(meta_s->data);
  Point *t = reinterpret_cast<Point *>(meta_t->data);
  Node *n = reinterpret_cast<Node *>(meta_g->data);

  int dim3 = pow(8, unif_level);

  // TODO: This bit where the message is interpreted might be made easier with
  // ReadBuffer
  int *meta = reinterpret_cast<int *>(static_cast<char *>(args));
  int first = meta[0];
  int last = meta[1];
  int range = last - first + 1;
  int recv_ns = meta[2];
  int recv_nt = meta[3];
  int *count_s = &meta[4]; // Used only if recv_ns > 0
  int *count_t = count_s + range * (recv_ns > 0); // Used only if recv_nt > 0
  Point *recv_s =
    reinterpret_cast<Point *>(static_cast<char *>(args) + sizeof(int) * 4 +
                              sizeof(int) * range * (recv_ns > 0) +
                              sizeof(int) * range * (recv_nt > 0));
  Point *recv_t = recv_s + recv_ns;

  hpx_addr_t done = hpx_lco_and_new(range * 2);

  if (recv_ns) {
    char type = 's';
    for (int i = first; i <= last; ++i) {
      Node *ns = &n[i];
      int incoming_ns = count_s[i - first];
      if (incoming_ns) {
        hpx_call(HPX_HERE, merge_points_action, done,
                 &s, &recv_s, &ns, &incoming_ns, &count[i], &type);
        recv_s += incoming_ns;
      } else {
        hpx_lco_and_set(done, HPX_NULL);
      }
    }
  } else {
    hpx_lco_and_set_num(done, range, HPX_NULL);
  }

  if (recv_nt) {
    char type = 't';
    for (int i = first; i <= last; ++i) {
      Node *nt = &n[i + dim3];
      int incoming_nt = count_t[i - first];
      if (incoming_nt) {
        hpx_call(HPX_HERE, merge_points_action, done,
                 &t, &recv_t, &nt, &incoming_nt, &count[i + dim3], &type);
        recv_t += incoming_nt;
      } else {
        hpx_lco_and_set(done, HPX_NULL);
      }
    }
  } else {
    hpx_lco_and_set_num(done, range, HPX_NULL);
  }

  // Wait until the data has been merged before releasing the parcel
  // NOTE: This 'done' will trigger once points are merged. The action that
  // triggers this will spawn more work, but not in a synchronous way.
  hpx_lco_wait(done);
  hpx_lco_delete_sync(done);
  hpx_lco_release(*user_lco, user_lco_buffer);
  hpx_gas_unpin(curr_unif_count);
  hpx_gas_unpin(curr_unif_count);
  hpx_gas_unpin(curr_sorted_src);
  hpx_gas_unpin(curr_sorted_tar);
  hpx_gas_unpin(curr_unif_grid);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED, recv_points_action,
           recv_points_handler, HPX_POINTER, HPX_SIZE_T);

// This is pretty straightforward. This rank will send to the given rank all
// those particles that are to be shipped out. There is nothing too complicated
// in this, just some indexing and so forth.
//
// NOTE: This would be a bit nicer looking using the Buffer types.
int send_points_handler(int rank, int *count_s, int *count_t,
                        int *offset_s, int *offset_t,
                        Point *sources, Point *targets) {
  // Note: all the pointers are local to the calling rank.
  int first = (rank == 0 ? 0 : distribute[rank - 1] + 1);
  int last = distribute[rank];
  int range = last - first + 1;
  int send_ns = 0, send_nt = 0;
  for (int i = first; i <= last; ++i) {
    send_ns += count_s[i];
    send_nt += count_t[i];
  }

  // Parcel message size
  size_t size = sizeof(int) * 4;
  if (send_ns) {
    size += sizeof(int) * range + sizeof(Point) * send_ns;
  }
  if (send_nt) {
    size += sizeof(int) * range + sizeof(Point) * send_nt;
  }

  // Acquire parcel
  hpx_parcel_t *p = hpx_parcel_acquire(NULL, size);
  void *data = hpx_parcel_get_data(p);
  int *meta = reinterpret_cast<int *>(static_cast<char *>(data));
  meta[0] = first;
  meta[1] = last;
  meta[2] = send_ns;
  meta[3] = send_nt;

  int *count = &meta[4];
  if (send_ns) {
    memcpy(count, &count_s[first], sizeof(int) * range);
    count += range;
  }

  if (send_nt) {
    memcpy(count, &count_t[first], sizeof(int) * range);
  }

  char *meta_s = static_cast<char *>(data) +
    sizeof(int) * (4 + range * (send_ns > 0) + range * (send_nt > 0));
  if (send_ns) {
    memcpy(meta_s, &sources[offset_s[first]], sizeof(Point) * send_ns);
  }

  char *meta_t = meta_s + send_ns * sizeof(Point);
  if (send_nt) {
    memcpy(meta_t, &targets[offset_t[first]], sizeof(Point) * send_nt);
  }

#if 0
  std::cout << send_ns << "\n" << std::flush;
#endif


  hpx_parcel_set_target(p, HPX_THERE(rank));
  hpx_parcel_set_action(p, recv_points_action);
  hpx_parcel_send(p, HPX_NULL);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, send_points_action, send_points_handler,
           HPX_INT, HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_POINTER,
           HPX_POINTER, HPX_POINTER);

int recv_node_handler(void *args, size_t size) {
  int *compressed_tree = reinterpret_cast<int *>(static_cast<char *>(args));
  int rank = hpx_get_my_rank();
  hpx_addr_t curr_unif_grid = hpx_addr_add(unif_grid,
                                           sizeof(ArrayMetaData) * rank,
                                           sizeof(ArrayMetaData));
  ArrayMetaData *meta;
  if (!hpx_gas_try_pin(curr_unif_grid, (void **)&meta))
    return HPX_ERROR;

  int type = compressed_tree[0];
  int id = compressed_tree[1];
  int n_nodes = compressed_tree[2];
  int dim3 = pow(8, unif_level);
  Node *n = reinterpret_cast<Node *>(meta->data);
  Node *curr = &n[id + type * dim3];

  if (n_nodes) {
    const int *branch = &compressed_tree[3];
    const int *tree = &compressed_tree[3 + n_nodes];
    curr->extract(branch, tree, n_nodes);
  }

  hpx_lco_and_set_num(curr->complete(), 8, HPX_NULL);

  hpx_gas_unpin(curr_unif_grid);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED, recv_node_action, recv_node_handler,
           HPX_POINTER, HPX_SIZE_T);

// This action is called once the
int send_node_handler(Node *n, ArrayMetaData *meta, int id, int type) {
  Node *curr = &n[id];
  int first = curr->first();
  int last = curr->last();
  int range = last - first + 1;
  Point *p = reinterpret_cast<Point *>(meta->data);

  // Rearrange points
  int *map = (type ? map_tar : map_src);
  Point *temp = new Point[range];
  for (int i = first; i <= last; ++i)
    temp[i - first] = p[map[i]];
  memcpy(p + first, temp, sizeof(Point) * range);
  delete [] temp;

  // Exclude curr as it is already allocated on remote localities
  int n_nodes = curr->n_descendants() - 1;
  int *compressed_tree = new int[3 + n_nodes * 2]();

  compressed_tree[0] = type; // source tree is 0, target tree is 1
  compressed_tree[1] = id; // where to merge
  compressed_tree[2] = n_nodes; // # of nodes
  if (n_nodes) {
    int *branch = &compressed_tree[3];
    int *tree = &compressed_tree[3 + n_nodes];
    int pos = 0;
    curr->compress(branch, tree, -1, pos);
  }

  int rank = hpx_get_my_rank();
  int num_ranks = hpx_get_num_ranks();
  hpx_addr_t done = hpx_lco_and_new(num_ranks - 1);

  for (int r = 0; r < num_ranks; ++r) {
    if (r != rank) {
      hpx_parcel_t *p = hpx_parcel_acquire(compressed_tree,
                                           sizeof(int) * (3 + n_nodes * 2));
      hpx_parcel_set_target(p, HPX_THERE(r));
      hpx_parcel_set_action(p, recv_node_action);
      hpx_parcel_send(p, done);
    }
  }

  // Clear local memory once all parcels are sent
  hpx_lco_wait(done);

  delete [] compressed_tree;
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, send_node_action, send_node_handler,
           HPX_POINTER, HPX_POINTER, HPX_INT, HPX_INT);


// This appears to be the main action that creates the trees. This will
// organize and call out to the other actions.
// NOTE: One thing is sure, this ought to be factored
int create_dual_tree_handler(hpx_addr_t sources, hpx_addr_t targets,
                             char exchange) {
  int rank = hpx_get_my_rank();
  int num_ranks = hpx_get_num_ranks();

  hpx_addr_t curr_s = hpx_addr_add(sources, sizeof(ArrayMetaData) * rank,
                                   sizeof(ArrayMetaData));
  hpx_addr_t curr_t = hpx_addr_add(targets, sizeof(ArrayMetaData) * rank,
                                   sizeof(ArrayMetaData));
  ArrayMetaData *meta_s{nullptr}, *meta_t{nullptr};
  if (!hpx_gas_try_pin(curr_s, (void **)&meta_s) ||
      !hpx_gas_try_pin(curr_t, (void **)&meta_t)) {
    return HPX_ERROR;
  }
  hpx_gas_unpin(curr_s);
  hpx_gas_unpin(curr_t);

  int n_sources = meta_s->count;
  int n_targets = meta_t->count;
  Point *p_s = reinterpret_cast<Point *>(meta_s->data);
  Point *p_t = reinterpret_cast<Point *>(meta_t->data);

  // Assign points to uniform grid
  int *gid_of_sources = new int[n_sources]();
  int *gid_of_targets = new int[n_targets]();
  int dim3 = pow(8, unif_level);
  int *local_count = new int[dim3 * 2]();
  int *local_scount = local_count;
  int *local_tcount = &local_count[dim3];
  // TODO: perhaps these are actions? That is a coarse parallelism - then
  // perhaps more might be added inside these functions
  assign_points_to_unif_grid(p_s, n_sources, gid_of_sources,
                             local_scount, 1.0 / size);
  assign_points_to_unif_grid(p_t, n_targets, gid_of_targets,
                             local_tcount, 1.0 / size);

  // Exchange counts
  //
  // TODO: Decide if it is better to just do a single reduction here. The
  // difference is between the current ~N^2 remote messages vs
  // ~2N remote messages. The only tradeoff is that the N^2 version might allow
  // lucky early ranks to start sooner.
  for (int r = 0; r < num_ranks; ++r) {
    hpx_parcel_t *p = hpx_parcel_acquire(local_count,
                                         sizeof(int) * dim3 * 2);
    hpx_parcel_set_target(p, HPX_THERE(r));
    hpx_parcel_set_action(p, exchange_count_action);
    hpx_parcel_send(p, HPX_NULL);
  }

  // Put points of the same grid together while waiting for
  // exchange_count_action to complete
  Point *temp_s = new Point[n_sources];
  Point *temp_t = new Point[n_targets];
  // TODO: Perhaps start these as actions to get at least that coarse
  // parallelism.
  int *local_offset_s = group_points_on_unif_grid(p_s, n_sources,
                                                  gid_of_sources, local_scount,
                                                  temp_s);
  int *local_offset_t = group_points_on_unif_grid(p_t, n_targets,
                                                  gid_of_targets, local_tcount,
                                                  temp_t);
  delete [] gid_of_sources;
  delete [] gid_of_targets;

  // Compute point distribution
  // First by getting this ranks LCO address
  hpx_addr_t curr_unif_count = hpx_addr_add(unif_count,
                                            sizeof(hpx_addr_t) * rank,
                                            sizeof(hpx_addr_t));
  hpx_addr_t *user_lco{nullptr};
  if (!hpx_gas_try_pin(curr_unif_count, (void **)&user_lco))
    return HPX_ERROR;

  // Second by getting a reference to the LCO data - if there is just the
  // one LCO, this would instead be a get probably.
  void *user_lco_buffer{nullptr};
  hpx_lco_getref(*user_lco, 1, &user_lco_buffer);
  int *global_count =
    reinterpret_cast<int *>(static_cast<char *>(user_lco_buffer)
                            + sizeof(int) * 2);

  // NOTE: inside this call is where the distribute array is allocated
  distribute = distribute_points(num_ranks, global_count, dim3);
  assert(distribute != nullptr);


  // Exchange points
  hpx_addr_t curr_sorted_s = hpx_addr_add(sorted_src,
                                          sizeof(ArrayMetaData) * rank,
                                          sizeof(ArrayMetaData));
  hpx_addr_t curr_sorted_t = hpx_addr_add(sorted_tar,
                                          sizeof(ArrayMetaData) * rank,
                                          sizeof(ArrayMetaData));
  hpx_addr_t curr_unif_grid = hpx_addr_add(unif_grid,
                                           sizeof(ArrayMetaData) * rank,
                                           sizeof(ArrayMetaData));
  hpx_addr_t curr_unif_done = hpx_addr_add(unif_done,
                                           sizeof(hpx_addr_t) * rank,
                                           sizeof(hpx_addr_t));

  ArrayMetaData *meta_g{nullptr};
  hpx_addr_t *gate{nullptr};

  if (!hpx_gas_try_pin(curr_sorted_s, (void **)&meta_s) ||
      !hpx_gas_try_pin(curr_sorted_t, (void **)&meta_t) ||
      !hpx_gas_try_pin(curr_unif_grid, (void **)&meta_g) ||
      !hpx_gas_try_pin(curr_unif_done, (void **)&gate))
    return HPX_ERROR;

  Node *ns = reinterpret_cast<Node *>(meta_g->data);
  Node *nt = reinterpret_cast<Node *>(meta_g->data + sizeof(Node) * dim3);

  // ?
  int *global_offset_s = init_point_exchange(rank, global_count, local_scount,
                                             local_offset_s, temp_s,
                                             ns, meta_s, 's');
  int *global_offset_t = init_point_exchange(rank, global_count + dim3,
                                             local_tcount, local_offset_t,
                                             temp_t, nt, meta_t, 't');
  hpx_lco_and_set(*gate, HPX_NULL);

  // So this one is pretty simple. It sends those points from this rank
  // going to the other rank in a parcel.
  for (int r = 0; r < num_ranks; ++r) {
    if (r != rank) {
      hpx_call(HPX_HERE, send_points_action, HPX_NULL, &r, &local_scount,
               &local_tcount, &local_offset_s, &local_offset_t,
               &temp_s, &temp_t);
    }
  }

  hpx_addr_t dual_tree_complete = hpx_lco_and_new(2 * dim3);

  if (exchange == 'y') {
    for (int r = 0; r < num_ranks; ++r) {
      int first = (r == 0 ? 0 : distribute[r - 1] + 1);
      int last = distribute[r];
      int s{0}, t{1};

      if (r == rank) {
        for (int i = first; i <= last; ++i) {
          if (global_count[i] == 0) {
            hpx_lco_and_set(dual_tree_complete, HPX_NULL);
          } else {
            hpx_call_when_with_continuation(ns[i].complete(), HPX_HERE,
                                            send_node_action, dual_tree_complete,
                                            hpx_lco_set_action, &ns, &meta_s,
                                            &i, &s);
          }

          if (global_count[i + dim3] == 0) {
            hpx_lco_and_set(dual_tree_complete, HPX_NULL);
          } else {
            hpx_call_when_with_continuation(nt[i].complete(), HPX_HERE,
                                            send_node_action, dual_tree_complete,
                                            hpx_lco_set_action, &nt, &meta_t,
                                            &i, &t);
          }
        }
      } else {
        for (int i = first; i <= last; ++i) {
          if (global_count[i] == 0) {
            hpx_lco_and_set(dual_tree_complete, HPX_NULL);
          } else {
            hpx_call_when(ns[i].complete(), dual_tree_complete,
                          hpx_lco_set_action, HPX_NULL, NULL, 0);
          }

          if (global_count[i + dim3] == 0) {
            hpx_lco_and_set(dual_tree_complete, HPX_NULL);
          } else {
            hpx_call_when(nt[i].complete(), dual_tree_complete,
                          hpx_lco_set_action, HPX_NULL, NULL, 0);
          }
        }
      }
    }
  } else {
    for (int r = 0; r < num_ranks; ++r) {
      int first = (r == 0 ? 0 : distribute[r - 1] + 1);
      int last = distribute[r];

      if (r == rank) {
        for (int i = first; i <= last; ++i) {
          if (global_count[i] == 0) {
            hpx_lco_and_set(dual_tree_complete, HPX_NULL);
          } else {
            hpx_call_when(ns[i].complete(), dual_tree_complete,
                          hpx_lco_set_action, HPX_NULL);
          }

          if (global_count[i + dim3] == 0) {
            hpx_lco_and_set(dual_tree_complete, HPX_NULL);
          } else {
            hpx_call_when(nt[i].complete(), dual_tree_complete,
                          hpx_lco_set_action, HPX_NULL);
          }
        }
      } else {
        for (int i = first; i <= last; ++i)
          hpx_lco_and_set_num(dual_tree_complete, 2, HPX_NULL);
      }
    }
  }

  hpx_lco_wait(dual_tree_complete);
  hpx_lco_delete_sync(dual_tree_complete);

  delete [] local_count;
  delete [] temp_s;
  delete [] temp_t;
  delete [] local_offset_s;
  delete [] local_offset_t;
  delete [] global_offset_s;
  delete [] global_offset_t;
  delete [] swap_src;
  delete [] bin_src;
  delete [] map_src;
  delete [] swap_tar;
  delete [] bin_tar;
  delete [] map_tar;

  hpx_lco_release(*user_lco, user_lco_buffer);
  hpx_gas_unpin(curr_unif_count);
  hpx_gas_unpin(curr_sorted_s);
  hpx_gas_unpin(curr_sorted_t);
  hpx_gas_unpin(curr_unif_grid);
  hpx_gas_unpin(curr_unif_done);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, create_dual_tree_action, create_dual_tree_handler,
           HPX_ADDR, HPX_ADDR, HPX_CHAR);

int finalize_partition_handler(void *unused, size_t size) {
  int rank = hpx_get_my_rank();
  hpx_addr_t count = hpx_addr_add(unif_count, sizeof(hpx_addr_t) * rank,
                                  sizeof(hpx_addr_t));
  hpx_addr_t done = hpx_addr_add(unif_done, sizeof(hpx_addr_t) * rank,
                                 sizeof(hpx_addr_t));
  hpx_addr_t grid = hpx_addr_add(unif_grid, sizeof(ArrayMetaData) * rank,
                                 sizeof(ArrayMetaData));
  hpx_addr_t sources = hpx_addr_add(sorted_src, sizeof(ArrayMetaData) * rank,
                                    sizeof(ArrayMetaData));
  hpx_addr_t targets = hpx_addr_add(sorted_tar, sizeof(ArrayMetaData) * rank,
                                    sizeof(ArrayMetaData));

  hpx_addr_t *user_lco{nullptr}, *gate{nullptr};
  ArrayMetaData *meta_g{nullptr}, *meta_s{nullptr}, *meta_t{nullptr};

  if (!hpx_gas_try_pin(count, (void **)&user_lco) ||
      !hpx_gas_try_pin(done, (void **)&gate) ||
      !hpx_gas_try_pin(grid, (void **)&meta_g) ||
      !hpx_gas_try_pin(sources, (void **)&meta_s) ||
      !hpx_gas_try_pin(targets, (void **)&meta_t)) {
    return HPX_ERROR;
  }

  hpx_lco_delete_sync(*user_lco);
  hpx_lco_delete_sync(*gate);
  delete [] meta_s->data;
  delete [] meta_t->data;

  int dim3 = pow(8, unif_level);
  Node *n = reinterpret_cast<Node *>(meta_g->data);

  int first = (rank == 0 ? 0 : distribute[rank - 1] + 1);
  int last = distribute[rank];

  for (int i = 0; i < first; ++i) {
    Node *curr = &n[i];
    hpx_lco_delete_sync(curr->sema());
    hpx_lco_delete_sync(curr->complete());

    for (int j = 0; j < 8; ++j) {
      Node *child = curr->child(j);
      if (child)
        child->destroy(true);
    }

    for (int j = 0; j < 8; ++j) {
      Node *child = curr->child(j);
      if (child) {
        delete [] child;
        break;
      }
    }
  }

  for (int i = first; i <= last; ++i) {
    Node *curr = &n[i];
    hpx_lco_delete_sync(curr->sema());
    hpx_lco_delete_sync(curr->complete());

    for (int j = 0; j < 8; ++j) {
      Node *child = curr->child(j);
      if (child)
        child->destroy(false);
    }
  }

  for (int i = last + 1; i < dim3; ++i) {
    Node *curr = &n[i];
    hpx_lco_delete_sync(curr->sema());
    hpx_lco_delete_sync(curr->complete());

    for (int j = 0; j < 8; ++j) {
      Node *child = curr->child(j);
      if (child)
        child->destroy(true);
    }

    for (int j = 0; j < 8; ++j) {
      Node *child = curr->child(j);
      if (child) {
        delete [] child;
        break;
      }
    }
  }

  delete [] meta_g->data;
  delete [] distribute;

  hpx_gas_unpin(count);
  hpx_gas_unpin(done);
  hpx_gas_unpin(grid);
  hpx_gas_unpin(sources);
  hpx_gas_unpin(targets);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED, finalize_partition_action,
           finalize_partition_handler, HPX_POINTER, HPX_SIZE_T);

/// The main action of the code starts here
int main_handler(char scaling, char datatype, char exchange,
                 int nsrc, int ntar, int threshold, int nseed) {
  int num_ranks = hpx_get_num_ranks();

  // Generate input data
  hpx_addr_t sources = hpx_gas_calloc_cyclic(num_ranks,
                                             sizeof(ArrayMetaData), 0);
  hpx_addr_t targets = hpx_gas_calloc_cyclic(num_ranks,
                                             sizeof(ArrayMetaData), 0);

  hpx_bcast_rsync(generate_input_action, &scaling, &datatype, &nsrc,
                  &ntar, &nseed, &sources, &targets);

  // Partition points to create dual tree

  hpx_time_t timer_start = hpx_time_now();

  // Determine domain geometry - this is likely better as a reduction LCO;
  //  the flexibility of the USER LCO is not needed.
  hpx_addr_t domain_geometry =
    hpx_lco_user_new(sizeof(int) + sizeof(double) * 6,
                     domain_geometry_init_action,
                     domain_geometry_op_action,
                     domain_geometry_predicate_action,
                     &num_ranks, sizeof(num_ranks));

  hpx_bcast_lsync(set_domain_geometry_action, HPX_NULL,
                  &sources, &targets, &domain_geometry);

  // NOTE: getref here is not a huge savings. It cuts down on the stack for
  // this thread a little. hpx_lco_get is likely easier
  void *temp{nullptr};
  hpx_lco_getref(domain_geometry, 1, &temp);
  double *var = reinterpret_cast<double *>(static_cast<char *>(temp) +
                                           sizeof(int));
  double size_x = var[1] - var[0];
  double size_y = var[3] - var[2];
  double size_z = var[5] - var[4];
  size = fmax(size_x, fmax(size_y, size_z));
  corner_x = (var[1] + var[0] - size) / 2;
  corner_y = (var[3] + var[2] - size) / 2;
  corner_z = (var[5] + var[4] - size) / 2;

  // Measure the reduction time
  hpx_time_t timer_domain = hpx_time_now();

  // Choose uniform partition level such that the number of grids is no less
  // than the number of ranks
  unif_level = ceil(log(num_ranks) / log(8)) + 1;

  unif_count = hpx_gas_calloc_cyclic(num_ranks, sizeof(hpx_addr_t), 0);
  unif_done = hpx_gas_calloc_cyclic(num_ranks, sizeof(hpx_addr_t), 0);
  unif_grid = hpx_gas_calloc_cyclic(num_ranks, sizeof(ArrayMetaData), 0);
  sorted_src = hpx_gas_calloc_cyclic(num_ranks, sizeof(ArrayMetaData), 0);
  sorted_tar = hpx_gas_calloc_cyclic(num_ranks, sizeof(ArrayMetaData), 0);

  hpx_bcast_rsync(init_partition_action, &unif_count, &unif_done, &unif_grid,
                  &sorted_src, &sorted_tar, &unif_level, &threshold,
                  &corner_x, &corner_y, &corner_z, &size);

  // Measure the partitioning setup. This is all just getting this and that
  // ready to do. No real work is done in the above.
  hpx_time_t timer_middle = hpx_time_now();

  hpx_bcast_rsync(create_dual_tree_action, &sources, &targets, &exchange);

  // All done, spit out some timing information before cleaning up and halting
  // TODO: improve the description once I read stuff
  hpx_time_t timer_end = hpx_time_now();
  double elapsed_total = hpx_time_diff_ms(timer_start, timer_end) / 1e3;
  double elapsed_reduction = hpx_time_diff_ms(timer_start, timer_domain) / 1e3;
  double elapsed_first = hpx_time_diff_ms(timer_domain, timer_middle) / 1e3;
  double elapsed_second = hpx_time_diff_ms(timer_middle, timer_end) / 1e3;
  std::cout << "Dual tree creation time: " << elapsed_total << "\n";
  std::cout << "  Reduction: " << elapsed_reduction << "\n";
  std::cout << "  Setup: " << elapsed_first << "\n";
  std::cout << "  Partition: " << elapsed_second << "\n";

  // This is just deleting the allocated resources. Nothing too special here.
  hpx_bcast_rsync(finalize_partition_action, NULL, 0);

  // Destroy input data
  hpx_bcast_rsync(destroy_input_action, &sources, &targets);

  hpx_lco_release(domain_geometry, temp);
  hpx_lco_delete_sync(domain_geometry);
  hpx_gas_free_sync(unif_count);
  hpx_gas_free_sync(unif_done);
  hpx_gas_free_sync(unif_grid);
  hpx_gas_free_sync(sorted_src);
  hpx_gas_free_sync(sorted_tar);
  hpx_gas_free_sync(sources);
  hpx_gas_free_sync(targets);

  hpx_exit(0, nullptr);
}
HPX_ACTION(HPX_DEFAULT, 0, main_action, main_handler,
           HPX_CHAR, HPX_CHAR, HPX_CHAR, HPX_INT, HPX_INT, HPX_INT, HPX_INT);

void Node::partition(Point *p, int *swap, int *bin, int *map,
                     int threshold, double corner_x, double corner_y,
                     double corner_z, double size) {
  int num_points = last_ - first_ + 1;
  assert(first_ >= 0 && last_ >= 0 && num_points >= 1);
  bool is_leaf = (num_points <= threshold ? true : false);

  if (parent_) {
    hpx_call_when(complete_, parent_->complete(), hpx_lco_set_action,
                  HPX_NULL, NULL, 0);
  }

  if (is_leaf) {
    // Set complete_ LCO of itself
    hpx_lco_and_set_num(complete_, 8, HPX_NULL);
  } else {
    // Compute center of the node
    double h = size / pow(2, idx_.level());
    double center_x = corner_x + (idx_.x() + 0.5) * h;
    double center_y = corner_y + (idx_.y() + 0.5) * h;
    double center_z = corner_z + (idx_.z() + 0.5) * h;

    // TODO --- this section should be factored probably
    //      --- though, the extra storage is baked into the interface
    // Sort points into 8 octants of the node
    int stat[8] = {0};   /// Count octant occupation
    int offset[8] = {0}; /// First position of each octant's points

    for (int i = first_; i <= last_; ++i) {
      Point *j = &p[map[i]];
      int temp = 4 * (j->z() > center_z) + 2 * (j->y() > center_y) +
        (j->x() > center_x);
      bin[i] = temp;
      stat[temp]++;
    }

    offset[0] = first_;
    for (int i = 1; i < 8; ++i) {
      offset[i] = offset[i - 1] + stat[i - 1];
    }

    // Rearrange map
    for (int i = first_; i <= last_; ++i) {
      int j = offset[bin[i]]++;
      swap[j] = map[i];
      // could just set bin[j] = map[i] here instead; the bin is no longer
      // used for the given i
    }

    for (int i = first_; i <= last_; ++i) {
      map[i] = swap[i];
    }
    // TODO --- end of the 'factor this' section

    // Create child nodes
    for (int i = 0; i < 8; ++i) {
      if (stat[i]) {
        Node *child = new Node{idx_.child(i), offset[i] - stat[i],
                               offset[i] - 1, this};
        child_[i] = child;

        hpx_call(HPX_HERE, partition_node_action, HPX_NULL,
                 &child, &p, &swap, &bin, &map);
      } else {
        hpx_lco_and_set(complete_, HPX_NULL);
      }
    }
  }
}

int Node::n_descendants() const {
  int count = 1;
  for (int i = 0; i < 8; ++i) {
    if (child_[i] != nullptr)
      count += child_[i]->n_descendants();
  }
  return count;
}

void Node::compress(int *branch, int *tree, int parent, int &curr) const {
  for (int i = 0; i < 8; ++i) {
    if (child_[i] != nullptr) {
      branch[curr] = i; // tracks which child exists
      tree[curr] = parent; // tracks the parent of the node being processed
      curr++; // Move onto the next slot
      // curr - 1 is the parent location for the subtree rooted at child_[i]
      child_[i]->compress(branch, tree, curr - 1, curr);
    }
  }
}

void Node::extract(const int *branch, const int *tree, int n_nodes) {
  // Extract a compressed remote tree representation. As the tree is remote,
  // only {parent_, child_, idx_} fields are needed.

  Node *descendants = new Node[n_nodes];

  // The compressed tree is created in depth first fashion. And there are two
  // choices here to fill in the parent_, child_, and idx_ fields of the
  // descendants.

  // Approach I: Setup parent_ and child_, which is an embarassingly parallel
  // operation on the @p branch and @p tree. Afterwards, fan out along the tree
  // to fill in idx.

  // Approach II: Go over the input sequentially. For each node encountered, by
  // the depth first property, the index of its parent is already set. So one
  // can finish in one loop.

  // If on each rank, there are multiple subtrees being merged, approach II
  // might be sufficient. Approach I can be considered if finer granularity is
  // needed.

  // Approach II is implemented here.
  for (int i = 0; i < n_nodes; ++i) {
    int pos = tree[i];
    int which = branch[i];
    Node *curr = &descendants[i];
    Node *parent = (pos < 0 ? this : &descendants[pos]);

    curr->set_parent(parent);
    curr->set_index(parent->index().child(which));
    parent->set_child(which, curr);
  }
}

void Node::destroy(bool allocated_in_array) {
  for (int i = 0; i < 8; ++i) {
    if (child_[i])
      child_[i]->destroy(allocated_in_array);
  }
  if (sema_ != HPX_NULL)
    hpx_lco_delete_sync(sema_);
  if (complete_ != HPX_NULL)
    hpx_lco_delete_sync(complete_);
  if (allocated_in_array == false)
    delete this;
}

