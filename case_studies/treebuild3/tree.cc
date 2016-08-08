#include "tree.h"

#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cassert>

#include <algorithm>
#include <functional>
#include <iostream>

#include "dashmm/domaingeometry.h"
#include "dashmm/reductionops.h"
#include "dashmm/shareddata.h"
using dashmm::int_sum_ident_op;
using dashmm::int_sum_op;
using dashmm::SharedData;
using dashmm::DomainGeometry;
using dashmm::Point;
using dashmm::Index;
using dashmm::Array;
using dashmm::ArrayRef;
using dashmm::ArrayData;


// The domain geometry
SharedData<DomainGeometry> domain{HPX_NULL};

int unif_level;         /// The level of uniform partition
hpx_addr_t unif_count;  /// The address of a reduction LCO used to perform the
                        /// counting of the source and target points
int *unif_count_value;  /// A local array holding the results
Node *unif_grid;        /// A rankwise global holding the uniform grid
hpx_addr_t unif_done;   /// A rankwise global holding an LCO

size_t sorted_src_count; /// How many sources this rank has
Point *sorted_src;       /// The sorted sources
size_t sorted_tar_count; /// How many targets this rank has
Point *sorted_tar;       /// The sorted targets
int *distribute;        /// Array giving the distribution of the uniform grid

int threshold;
int *swap_src;
int *bin_src;
int *map_src;
int *swap_tar;
int *bin_tar;
int *map_tar;

/////////////////////////////////////////////////////////////////////
// A couple general routines
/////////////////////////////////////////////////////////////////////

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


/////////////////////////////////////////////////////////////////////
// Domain reduction stuff
/////////////////////////////////////////////////////////////////////

int set_domain_geometry_handler(hpx_addr_t sources_gas,
                                hpx_addr_t targets_gas,
                                hpx_addr_t dom_gas,
                                hpx_addr_t domain_geometry) {
  domain = SharedData<DomainGeometry>{dom_gas};

  Array<Point> sources{sources_gas};
  ArrayRef<Point> src_ref = sources.ref();
  ArrayData<Point> src_data = src_ref.pin();
  Point *s = src_data.value();

  Array<Point> targets{targets_gas};
  ArrayRef<Point> trg_ref = targets.ref();
  ArrayData<Point> trg_data = trg_ref.pin();
  Point *t = trg_data.value();

  double var[6] = {1e50, -1e50, 1e50, -1e50, 1e50, -1e50};

  // NOTE: Here is an opportunity to do even more parallelism. This will be
  // a single thread per locality. Why not do the on locality reduction in
  // parallel.
  for (size_t i = 0; i < src_ref.n(); ++i) {
    var[0] = fmin(var[0], s[i].x());
    var[1] = fmax(var[1], s[i].x());
    var[2] = fmin(var[2], s[i].y());
    var[3] = fmax(var[3], s[i].y());
    var[4] = fmin(var[4], s[i].z());
    var[5] = fmax(var[5], s[i].z());
  }

  for (size_t i = 0; i < trg_ref.n(); ++i) {
    var[0] = fmin(var[0], t[i].x());
    var[1] = fmax(var[1], t[i].x());
    var[2] = fmin(var[2], t[i].y());
    var[3] = fmax(var[3], t[i].y());
    var[4] = fmin(var[4], t[i].z());
    var[5] = fmax(var[5], t[i].z());
  }

  hpx_lco_set_lsync(domain_geometry, sizeof(double) * 6, var, HPX_NULL);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, set_domain_geometry_action,
           set_domain_geometry_handler,
           HPX_ADDR, HPX_ADDR, HPX_ADDR, HPX_ADDR);

void domain_geometry_init_handler(double *values, const size_t UNUSED) {
  values[0] = 1e50; // xmin
  values[1] = -1e50; // xmax
  values[2] = 1e50; // ymin
  values[3] = -1e50; // ymax
  values[4] = 1e50; // zmin
  values[5] = -1e50; // zmax
}
HPX_ACTION(HPX_FUNCTION, 0, domain_geometry_init_action,
           domain_geometry_init_handler, HPX_POINTER, HPX_SIZE_T);

void domain_geometry_op_handler(double *lhs, double *rhs, size_t UNUSED) {
  lhs[0] = fmin(lhs[0], rhs[0]);
  lhs[1] = fmax(lhs[1], rhs[1]);
  lhs[2] = fmin(lhs[2], rhs[2]);
  lhs[3] = fmax(lhs[3], rhs[3]);
  lhs[4] = fmin(lhs[4], rhs[4]);
  lhs[5] = fmax(lhs[5], rhs[5]);
}
HPX_ACTION(HPX_FUNCTION, 0, domain_geometry_op_action,
           domain_geometry_op_handler, HPX_POINTER, HPX_POINTER, HPX_SIZE_T);

void compute_domain_geometry(Array<Point> sources, Array<Point> targets) {
  // Allocate the shared data
  domain = SharedData<DomainGeometry>{nullptr};

  // Create a reduction LCO
  hpx_addr_t domain_geometry =
    hpx_lco_reduce_new(hpx_get_num_ranks(), sizeof(double) * 6,
                       domain_geometry_init_action,
                       domain_geometry_op_action);

  // Launch the reduction actions
  hpx_addr_t sglob = sources.data();
  hpx_addr_t tglob = targets.data();
  hpx_addr_t dglob = domain.data();
  hpx_bcast_lsync(set_domain_geometry_action, HPX_NULL,
                  &sglob, &tglob, &dglob, &domain_geometry);

  // Get the result
  double var[6];
  hpx_lco_get(domain_geometry, sizeof(double) * 6, &var);
  hpx_lco_delete_sync(domain_geometry);

  // Setup DomainGeometry
  double length = fmax(var[1] - var[0],
                       fmax(var[3] - var[2], var[5] - var[4]));
  DomainGeometry geo{Point{(var[1] + var[0] - length) / 2,
                           (var[3] + var[2] - length) / 2,
                           (var[5] + var[4] - length) / 2}, length};

  // Share with everyone
  domain.reset(&geo);
}

/////////////////////////////////////////////////////////////////////
// Basic Setup Stuff
/////////////////////////////////////////////////////////////////////

int init_partition_handler(hpx_addr_t count, int level, int limit) {
  int dim = pow(2, level), dim3 = pow(8, level);

  unif_count = count;
  unif_level = level;
  threshold = limit;

  // We here allocate space for the result of the counting
  unif_count_value = new int[dim3 * 2]();

  // Setup unif_done LCO
  unif_done = hpx_lco_and_new(1);

  // Setup unif_grid
  unif_grid = new Node[2 * dim3];
  for (int iz = 0; iz < dim; ++iz) {
    for (int iy = 0; iy < dim; ++iy) {
      for (int ix = 0; ix < dim; ++ix) {
        uint64_t mid = morton_key(ix, iy, iz);
        unif_grid[mid] = Node{Index{ix, iy, iz, level}};
        unif_grid[mid + dim3] = Node{Index{ix, iy, iz, level}};
      }
    }
  }

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, init_partition_action, init_partition_handler,
           HPX_ADDR, HPX_INT, HPX_INT);

void setup_basic_data(int thresh) {
  int num_ranks = hpx_get_num_ranks();

  // Choose uniform partition level such that the number of grids is no less
  // than the number of ranks
  int level = ceil(log(num_ranks) / log(8)) + 1;
  int dim3 = pow(8, unif_level);
  hpx_addr_t ucount = hpx_lco_reduce_new(num_ranks, sizeof(int) * (dim3 * 2),
                                        int_sum_ident_op,
                                        int_sum_op);

  hpx_bcast_rsync(init_partition_action, &ucount, &level, &thresh);
}

/////////////////////////////////////////////////////////////////////
// The tree creation work
/////////////////////////////////////////////////////////////////////

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

// This will assign the points to the uniform grid. This gives the points
// the id (in the Morton Key sense) of the box to which they are assigned,
// and it will count the numbers in each box.
void assign_points_to_unif_grid(const Point *P, int npts, int *gid,
                                int *count, const DomainGeometry *geo) {
  Point corner = geo->low();
  double scale = 1.0 / geo->size();

  // TODO: This is serial processing; is there some way to parallelize this?
  //   This would perhaps be worth timing.
  int dim = pow(2, unif_level);
  for (int i = 0; i < npts; ++i) {
    const Point *p = &P[i];
    int xid = std::min(dim - 1, (int)(dim * (p->x() - corner.x()) * scale));
    int yid = std::min(dim - 1, (int)(dim * (p->y() - corner.y()) * scale));
    int zid = std::min(dim - 1, (int)(dim * (p->z() - corner.z()) * scale));
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
  auto geo = domain.value();
  n->partition(p, swap, bin, map, threshold, geo->low(), geo->size());
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
                         const Point *temp, Node *n,
                         size_t sorted_count, Point **sorted,
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

  sorted_count = num_points;

  if (num_points > 0) {
    *sorted = new Point[num_points]();

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
        memcpy(*sorted + global_offset[i - first], &temp[local_offset[i]],
               sizeof(Point) * local_count[i]);

        if (local_count[i] == global_count[i]) {
          // This grid does not expect remote points.
          // Spawn adaptive partitioning
          hpx_call(HPX_HERE, partition_node_action, HPX_NULL,
                   &curr, sorted, swap, bin, map);
        } else {
          // Now use first_ as an iterator to trace the position to merge
          // remote points
          curr->set_first(global_offset[i - first] + local_count[i]);
        }
      }
    }
  } else {
    *sorted = nullptr;
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
int recv_points_handler(void *args, size_t UNUSED) {
  // Wait until the buffer is allocated before merging incoming messages
  // We could do this as a call when, but then we need to be aware of the
  // addresses for every rank's LCO. For now, we do this, as it is simpler.
  hpx_lco_wait(unif_done);

  // Now merge incoming message
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
      Node *ns = &unif_grid[i];
      int incoming_ns = count_s[i - first];
      if (incoming_ns) {
        hpx_call(HPX_HERE, merge_points_action, done,
                 &sorted_src, &recv_s, &ns, &incoming_ns,
                 &unif_count_value[i], &type);
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
      Node *nt = &unif_grid[i + dim3];
      int incoming_nt = count_t[i - first];
      if (incoming_nt) {
        hpx_call(HPX_HERE, merge_points_action, done,
                 &sorted_tar, &recv_t, &nt, &incoming_nt,
                 &unif_count_value[i + dim3], &type);
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

  // Parcel message length
  size_t bytes = sizeof(int) * 4;
  if (send_ns) {
    bytes += sizeof(int) * range + sizeof(Point) * send_ns;
  }
  if (send_nt) {
    bytes += sizeof(int) * range + sizeof(Point) * send_nt;
  }

  // Acquire parcel
  hpx_parcel_t *p = hpx_parcel_acquire(nullptr, bytes);
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

  hpx_parcel_set_target(p, HPX_THERE(rank));
  hpx_parcel_set_action(p, recv_points_action);
  hpx_parcel_send(p, HPX_NULL);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, send_points_action, send_points_handler,
           HPX_INT, HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_POINTER,
           HPX_POINTER, HPX_POINTER);

// This is the action on the other side that receives the partitioned tree.
int recv_node_handler(int *compressed_tree, size_t UNUSED) {
  int type = compressed_tree[0];
  int id = compressed_tree[1];
  int n_nodes = compressed_tree[2];
  int dim3 = pow(8, unif_level);
  Node *curr = &unif_grid[id + type * dim3];

  if (n_nodes) {
    const int *branch = &compressed_tree[3];
    const int *tree = &compressed_tree[3 + n_nodes];
    curr->extract(branch, tree, n_nodes);
  }

  hpx_lco_and_set_num(curr->complete(), 8, HPX_NULL);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED, recv_node_action, recv_node_handler,
           HPX_POINTER, HPX_SIZE_T);

// This action is called once the individual grids are done.
// Also, this is where the points are finally rearranged. All other work has
// been to set up their eventual index in the sorted situation. Here is the
// actual shuffling.
//
// This will send one message for each grid box to all other localities.
// This does allow for maximum parallelism. It is likely the right approach.
int send_node_handler(Node *n, Point *sorted, int id, int type) {
  Node *curr = &n[id];
  int first = curr->first();
  int last = curr->last();
  int range = last - first + 1;

  // Rearrange points
  int *map = (type ? map_tar : map_src);
  Point *temp = new Point[range];
  for (int i = first; i <= last; ++i) {
    temp[i - first] = sorted[map[i]];
  }
  memcpy(sorted + first, temp, sizeof(Point) * range);
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
//
// TODO: You are working through this. We have to stop briefly to see which
// of the other Arrays are actually needed in the GAS.
//
int create_dual_tree_handler(hpx_addr_t sources_gas, hpx_addr_t targets_gas) {
  int rank = hpx_get_my_rank();
  int num_ranks = hpx_get_num_ranks();

  Array<Point> sources{sources_gas};
  ArrayRef<Point> src_ref = sources.ref(rank);
  ArrayData<Point> src_data = src_ref.pin();
  Point *p_s = src_data.value();
  int n_sources = src_ref.n();

  Array<Point> targets{targets_gas};
  ArrayRef<Point> trg_ref = targets.ref(rank);
  ArrayData<Point> trg_data = trg_ref.pin();
  Point *p_t = trg_data.value();
  int n_targets = trg_ref.n();

  // Assign points to uniform grid
  int *gid_of_sources = new int[n_sources]();
  int *gid_of_targets = new int[n_targets]();
  int dim3 = pow(8, unif_level);
  int *local_count = new int[dim3 * 2]();
  int *local_scount = local_count;
  int *local_tcount = &local_count[dim3];
  // TODO: perhaps these are actions? That is a coarse parallelism - then
  // perhaps more might be added inside these functions
  auto geo = domain.value();
  assign_points_to_unif_grid(p_s, n_sources, gid_of_sources,
                             local_scount, geo.value());
  assign_points_to_unif_grid(p_t, n_targets, gid_of_targets,
                             local_tcount, geo.value());

  // Exchange counts
  hpx_lco_set(unif_count, sizeof(int) * dim3 * 2, local_count,
              HPX_NULL, HPX_NULL);

  // Put points of the same grid together while waiting for
  // counting to complete
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
  hpx_lco_get(unif_count, sizeof(int) * (dim3 * 2), unif_count_value);
  distribute = distribute_points(num_ranks, unif_count_value, dim3);
  assert(distribute != nullptr);


  // Exchange points
  Node *ns = unif_grid;
  Node *nt = &unif_grid[dim3];

  int *global_offset_s = init_point_exchange(rank, unif_count_value,
                                             local_scount, local_offset_s,
                                             temp_s, ns, sorted_src_count,
                                             &sorted_src, 's');
  int *global_offset_t = init_point_exchange(rank, unif_count_value + dim3,
                                             local_tcount, local_offset_t,
                                             temp_t, nt, sorted_tar_count,
                                             &sorted_tar, 't');
  hpx_lco_and_set(unif_done, HPX_NULL);

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

  for (int r = 0; r < num_ranks; ++r) {
    int first = (r == 0 ? 0 : distribute[r - 1] + 1);
    int last = distribute[r];
    int s{0}, t{1};

    if (r == rank) {
      for (int i = first; i <= last; ++i) {
        if (unif_count_value[i] == 0) {
          hpx_lco_and_set(dual_tree_complete, HPX_NULL);
        } else {
          hpx_call_when_with_continuation(ns[i].complete(), HPX_HERE,
                                          send_node_action, dual_tree_complete,
                                          hpx_lco_set_action, &ns, &sorted_src,
                                          &i, &s);
        }

        if (unif_count_value[i + dim3] == 0) {
          hpx_lco_and_set(dual_tree_complete, HPX_NULL);
        } else {
          hpx_call_when_with_continuation(nt[i].complete(), HPX_HERE,
                                          send_node_action, dual_tree_complete,
                                          hpx_lco_set_action, &nt, &sorted_tar,
                                          &i, &t);
        }
      }
    } else {
      for (int i = first; i <= last; ++i) {
        if (unif_count_value[i] == 0) {
          hpx_lco_and_set(dual_tree_complete, HPX_NULL);
        } else {
          hpx_call_when(ns[i].complete(), dual_tree_complete,
                        hpx_lco_set_action, HPX_NULL, NULL, 0);
        }

        if (unif_count_value[i + dim3] == 0) {
          hpx_lco_and_set(dual_tree_complete, HPX_NULL);
        } else {
          hpx_call_when(nt[i].complete(), dual_tree_complete,
                        hpx_lco_set_action, HPX_NULL, NULL, 0);
        }
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

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, create_dual_tree_action, create_dual_tree_handler,
           HPX_ADDR, HPX_ADDR);

void create_dual_tree(Array<Point> &sources, Array<Point> &targets) {
  hpx_addr_t source_gas = sources.data();
  hpx_addr_t target_gas = targets.data();
  hpx_bcast_rsync(create_dual_tree_action, &source_gas, &target_gas);
}

/////////////////////////////////////////////////////////////////////
// Finalize the partition work - that is, clean up
/////////////////////////////////////////////////////////////////////

int finalize_partition_handler(void *unused, size_t UNUSED) {
  int rank = hpx_get_my_rank();

  delete [] unif_count_value;
  hpx_lco_delete_sync(unif_done);
  delete [] sorted_src;
  delete [] sorted_tar;

  int dim3 = pow(8, unif_level);

  int first = (rank == 0 ? 0 : distribute[rank - 1] + 1);
  int last = distribute[rank];

  // TODO: Wait! Don't we need to delete the rest of the nodes too?
  // this should have had 2 * dim3 nodes in it.

  // NOTE: The difference here is that there are two different allocation
  // schemes for the nodes.

  for (int i = 0; i < first; ++i) {
    Node *curr = &unif_grid[i];
    hpx_lco_delete_sync(curr->sema());
    hpx_lco_delete_sync(curr->complete());

    for (int j = 0; j < 8; ++j) {
      Node *child = curr->child(j);
      if (child) {
        child->destroy(true);
      }
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
    Node *curr = &unif_grid[i];
    hpx_lco_delete_sync(curr->sema());
    hpx_lco_delete_sync(curr->complete());

    for (int j = 0; j < 8; ++j) {
      Node *child = curr->child(j);
      if (child) {
        child->destroy(false);
      }
    }
  }

  for (int i = last + 1; i < dim3; ++i) {
    Node *curr = &unif_grid[i];
    hpx_lco_delete_sync(curr->sema());
    hpx_lco_delete_sync(curr->complete());

    for (int j = 0; j < 8; ++j) {
      Node *child = curr->child(j);
      if (child) {
        child->destroy(true);
      }
    }

    for (int j = 0; j < 8; ++j) {
      Node *child = curr->child(j);
      if (child) {
        delete [] child;
        break;
      }
    }
  }

  delete [] unif_grid;
  delete [] distribute;

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED, finalize_partition_action,
           finalize_partition_handler, HPX_POINTER, HPX_SIZE_T);

void finalize_partition() {
  hpx_bcast_rsync(finalize_partition_action, NULL, 0);
  hpx_lco_delete_sync(unif_count);
}

/////////////////////////////////////////////////////////////////////

void Node::partition(Point *p, int *swap, int *bin, int *map,
                     int threshold, Point corner, double size) {
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
    double center_x = corner.x() + (idx_.x() + 0.5) * h;
    double center_y = corner.y() + (idx_.y() + 0.5) * h;
    double center_z = corner.z() + (idx_.z() + 0.5) * h;

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

