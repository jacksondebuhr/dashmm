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
                                hpx_addr_t domain_geometry) {
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
  // parallel. Likely the decision to do this will depend somewhat on how
  // many threads and how many points we have. Some smarts about this would
  // be nice to build in. However, the main case of interest is when there
  // are loads of particles per rank, so perhaps doing in a bunch of chunks
  // makes sense as a default, and use a lower bound as a cutoff or something.
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
           HPX_ADDR, HPX_ADDR, HPX_ADDR);

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

/// Compute the bounding box for the given source and target points.
///
/// This will return an LCO which will contain six doubles, the first three
/// being the low corner, and the second three being the high corner. This
/// is an asynchronous call, it will return as soon as the work is scheduled.
///
/// \param sources - the source array
/// \param targets - the target array
///
/// \returns - address of an LCO containing the reduced domain
hpx_addr_t compute_domain_geometry(Array<Point> sources,
                                   Array<Point> targets) {
  // Create a reduction LCO
  hpx_addr_t domain_geometry =
    hpx_lco_reduce_new(hpx_get_num_ranks(), sizeof(double) * 6,
                       domain_geometry_init_action,
                       domain_geometry_op_action);

  // Launch the reduction actions
  hpx_addr_t sglob = sources.data();
  hpx_addr_t tglob = targets.data();
  hpx_bcast_lsync(set_domain_geometry_action, HPX_NULL,
                  &sglob, &tglob, &domain_geometry);

  return domain_geometry;
}

/////////////////////////////////////////////////////////////////////
// Basic Setup Stuff
/////////////////////////////////////////////////////////////////////

int init_partition_handler(hpx_addr_t rwdata, hpx_addr_t count, int limit,
                           hpx_addr_t domain_geometry) {
  RankWise<DualTree> global_tree{rwdata};
  auto tree = global_tree.here();

  int num_ranks = hpx_get_num_ranks();
  tree->set_unif_level(ceil(log(num_ranks) / log(8)) + 1);
  int dim = pow(2, tree->unif_level());
  tree->set_dim3(pow(8, tree->unif_level()));
  tree->set_unif_count(count);
  tree->set_threshold(limit);

  // We here allocate space for the result of the counting
  tree->set_unif_count_value(new int[tree->dim3() * 2]());

  // Setup unif_done LCO
  tree->set_unif_done(hpx_lco_and_new(1));

  // Setup unif_grid
  Node *unif_grid = new Node[2 * tree->dim3()];
  for (int iz = 0; iz < dim; ++iz) {
    for (int iy = 0; iy < dim; ++iy) {
      for (int ix = 0; ix < dim; ++ix) {
        uint64_t mid = morton_key(ix, iy, iz);
        unif_grid[mid] = Node{Index{ix, iy, iz, tree->unif_level()}};
        unif_grid[mid + tree->dim3()] =
            Node{Index{ix, iy, iz, tree->unif_level()}};
      }
    }
  }
  tree->set_unif_grid(unif_grid);

  // Setup domain_
  double var[6];
  hpx_lco_get(domain_geometry, sizeof(double) * 6, &var);
  double length = fmax(var[1] - var[0],
                       fmax(var[3] - var[2], var[5] - var[4]));
  DomainGeometry geo{Point{(var[1] + var[0] - length) / 2,
                           (var[3] + var[2] - length) / 2,
                           (var[5] + var[4] - length) / 2}, length};
  tree->set_domain(geo);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, init_partition_action, init_partition_handler,
           HPX_ADDR, HPX_ADDR, HPX_INT, HPX_ADDR);

RankWise<DualTree> setup_basic_data(int threshold,
                                    hpx_addr_t domain_geometry) {
  RankWise<DualTree> retval{};
  retval.allocate();
  if (!retval.valid()) {
    // We return the invalid value to indicate the error
    return retval;
  }

  // Now the single things are created.
  int num_ranks = hpx_get_num_ranks();
  int level = ceil(log(num_ranks) / log(8)) + 1;
  int dim3 = pow(8, level);
  hpx_addr_t ucount = hpx_lco_reduce_new(num_ranks, sizeof(int) * (dim3 * 2),
                                        int_sum_ident_op,
                                        int_sum_op);
  hpx_addr_t rwdata = retval.data();
  hpx_bcast_rsync(init_partition_action, &rwdata, &ucount, &threshold,
                  &domain_geometry);

  return retval;
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
void assign_points_to_unif_grid(const Point *P, int npts,
                                const DomainGeometry &geo, int unif_level,
                                int *gid, int *count) {
  Point corner = geo.low();
  double scale = 1.0 / geo.size();

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
// TODO: Note that this too is a serial operation. Could we do some on-rank
// parallelism here?
//
// TODO: Work out what the correct way to parallelize this would be. And then
// make that happen.
int *group_points_on_unif_grid(Point *p_in, int npts, int dim3,
                               int *gid_of_points, const int *count) {
  int *offset = new int[dim3]();

  offset[0] = 0;
  for (int i = 1; i < dim3; ++i) {
    offset[i] = offset[i - 1] + count[i - 1];
  }

  // Set gid to the final location of the particle; NOTE: this will modify
  // the offset. This is corrected below
  for (int i = 0; i < npts; ++i) {
    int gid = gid_of_points[i];
    gid_of_points[i] = offset[gid];
    offset[gid]++;
  }

  // Declare some loop variables
  Point source, save;
  int source_sort, save_sort;
  int isource, isave, idest;

  // Do an O(N) rearrangement
  for (int i = 0; i < npts; ++i) {
    if (gid_of_points[i] != i) {
      source = p_in[i];
      source_sort = gid_of_points[i];
      isource = gid_of_points[i];
      idest = gid_of_points[i];

      do {
        save = p_in[idest];
        save_sort = gid_of_points[idest];
        isave = gid_of_points[idest];

        p_in[idest] = source;
        gid_of_points[idest] = source_sort;

        if (idest == i) break;

        source = save;
        source_sort = save_sort;
        isource = isave;

        idest = isource;
      } while (1);
    }
  }

  // Correct offset
  for (int i = 0; i < dim3; ++i) {
    offset[i] -= count[i];
  }

  return offset;
}


// New factor
int *sort_local_points(DualTree *tree, Point *p_s,
                       int n_sources, Point *p_t, int n_targets,
                       int **local_offset_s, int **local_offset_t) {
  int *local_count = new int[tree->dim3() * 2]();
  int *local_scount = local_count;
  int *local_tcount = &local_count[tree->dim3()];

  int *gid_of_sources = new int[n_sources]();
  int *gid_of_targets = new int[n_targets]();
  // TODO: perhaps these are actions? That is a coarse parallelism - then
  // perhaps more might be added inside these functions
  assign_points_to_unif_grid(p_s, n_sources, tree->domain(),
                             tree->unif_level(), gid_of_sources, local_scount);
  assign_points_to_unif_grid(p_t, n_targets, tree->domain(),
                             tree->unif_level(), gid_of_targets, local_tcount);

  // Exchange counts
  hpx_lco_set(tree->unif_count(), sizeof(int) * tree->dim3() * 2, local_count,
              HPX_NULL, HPX_NULL);

  // Put points of the same grid together while waiting for
  // counting to complete
  // TODO: Perhaps start these as actions to get at least that coarse
  // parallelism.
  *local_offset_s = group_points_on_unif_grid(p_s, n_sources, tree->dim3(),
                                                  gid_of_sources, local_scount);
  *local_offset_t = group_points_on_unif_grid(p_t, n_targets, tree->dim3(),
                                                  gid_of_targets, local_tcount);
  delete [] gid_of_sources;
  delete [] gid_of_targets;

  return local_count;
}


/// This is a thin wrapper around the partition method of a node. This is
/// made into an action to make use of the parallelism available with HPX.
int partition_node_handler(Node *n, Point *p, DomainGeometry *geo,
                           int threshold) {
  n->partition(p, threshold, geo);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, partition_node_action, partition_node_handler,
           HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_INT);

// This sets up some orgnaizational structures. The return value is an
// array of offsets in the final set of points for each part of the uniform
// grid. This is basically just setup. There is a chance that some work occurs
// in one branch. Before leaving, this will bring the points that do not have
// to change rank into their correct location. If it is detected that all of
// the points for that part of the tree have arrived, then the partitioning
// work will begin.
int *init_point_exchange(int rank, DualTree *tree,
                         const int *local_count, const int *local_offset,
                         const Point *temp, Node *n,
                         char type) {
  int first = tree->first(rank);
  int last = tree->last(rank);
  int range = last - first + 1;

  Point *sorted{nullptr};
  int *global_count{nullptr};

  if (type == 's') {
    global_count = tree->unif_count_value();
  } else {
    global_count = tree->unif_count_value() + tree->dim3();
  }

  // Compute global_offset
  int *global_offset = new int[range]();
  int num_points = global_count[first];
  for (int i = first + 1; i <= last; ++i) {
    num_points += global_count[i];
    global_offset[i - first] = global_offset[i - first - 1] +
                               global_count[i - 1];
  }
  size_t sorted_count = num_points;

  if (num_points > 0) {
    sorted = new Point[num_points]();

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
        memcpy(sorted + global_offset[i - first], &temp[local_offset[i]],
               sizeof(Point) * local_count[i]);

        if (local_count[i] == global_count[i]) {
          // This grid does not expect remote points.
          // Spawn adaptive partitioning
          // TODO: This is a mild cheat for the type...
          const DomainGeometry *arg = &(tree->domain());
          int threshold = tree->threshold();
          hpx_call(HPX_HERE, partition_node_action, HPX_NULL,
                   &curr, &sorted, &arg, &threshold);
        } else {
          // Now use first_ as an iterator to trace the position to merge
          // remote points
          curr->set_first(global_offset[i - first] + local_count[i]);
        }
      }
    }
  } else {
    sorted = nullptr;
  }

  if (type == 's') {
    tree->set_sorted_src_count(sorted_count);
    tree->set_sorted_src(sorted);
  } else {
    tree->set_sorted_tar_count(sorted_count);
    tree->set_sorted_tar(sorted);
  }

  return global_offset;
}

// This action merges particular points with the sorted list. Also, if this
// is the last set of points that are merged, this will go ahead and start
// the adaptive partitioning of that part of the local tree.
int merge_points_handler(Point *p, Point *temp, Node *n,
                         int n_arrived, int n_total, char type,
                         hpx_addr_t rwgas) {
  RankWise<DualTree> global_tree{rwgas};
  auto local_tree = global_tree.here();

  // Note: all the pointers are local to the calling rank.
  hpx_lco_sema_p(n->sema());
  int first = n->first();
  memcpy(p + first, temp, sizeof(Point) * n_arrived);
  n->set_first(first + n_arrived);

  if (first + n_arrived > n->last()) {
    // Spawn adaptive partitioning
    n->set_first(n->last() + 1 - n_total);
    // TODO: This is mildly a cheat...
    const DomainGeometry *geoarg = &(local_tree->domain());
    int thresh = local_tree->threshold();

    if (type == 's') {
      hpx_call(HPX_HERE, partition_node_action, HPX_NULL,
               &n, &p, &geoarg, &thresh);
    } else {
      hpx_call(HPX_HERE, partition_node_action, HPX_NULL,
               &n, &p, &geoarg, &thresh);
    }
  }
  hpx_lco_sema_v(n->sema(), HPX_NULL);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, merge_points_action, merge_points_handler,
           HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_INT, HPX_INT, HPX_CHAR,
           HPX_ADDR);

// This is the 'far-side' of the send points message. This action merges the
// incoming points into the sorted list and will then spawn the adaptive
// partition if this happens to be the last block for a given uniform grid.
int recv_points_handler(void *args, size_t UNUSED) {
  hpx_addr_t *rwarg = static_cast<hpx_addr_t *>(args);
  RankWise<DualTree> global_tree{*rwarg};
  auto local_tree = global_tree.here();

  // Wait until the buffer is allocated before merging incoming messages
  // We could do this as a call when, but then we need to be aware of the
  // addresses for every rank's LCO. For now, we do this, as it is simpler.
  hpx_lco_wait(local_tree->unif_done());

  // TODO: This bit where the message is interpreted might be made easier with
  // ReadBuffer
  int *meta = reinterpret_cast<int *>(static_cast<char *>(args)
                                      + sizeof(hpx_addr_t));
  int first = meta[0];
  int last = meta[1];
  int range = last - first + 1;
  int recv_ns = meta[2];
  int recv_nt = meta[3];
  int *count_s = &meta[4]; // Used only if recv_ns > 0
  int *count_t = count_s + range * (recv_ns > 0); // Used only if recv_nt > 0
  Point *recv_s =
    reinterpret_cast<Point *>(static_cast<char *>(args) + sizeof(int) * 4 +
                              sizeof(hpx_addr_t) +
                              sizeof(int) * range * (recv_ns > 0) +
                              sizeof(int) * range * (recv_nt > 0));
  Point *recv_t = recv_s + recv_ns;

  hpx_addr_t done = hpx_lco_and_new(range * 2);

  if (recv_ns) {
    char type = 's';
    for (int i = first; i <= last; ++i) {
      Node *ns = &(local_tree->unif_grid()[i]);
      int incoming_ns = count_s[i - first];
      if (incoming_ns) {
        Point *ss = local_tree->sorted_src();
        hpx_call(HPX_HERE, merge_points_action, done,
                 &ss, &recv_s, &ns, &incoming_ns,
                 &(local_tree->unif_count_value()[i]), &type, rwarg);
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
      Node *nt = &(local_tree->unif_grid()[i + local_tree->dim3()]);
      int incoming_nt = count_t[i - first];
      if (incoming_nt) {
        Point *st = local_tree->sorted_tar();
        hpx_call(HPX_HERE, merge_points_action, done,
                 &st, &recv_t, &nt, &incoming_nt,
                 &(local_tree->unif_count_value()[i + local_tree->dim3()]),
                 &type, rwarg);
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
                        Point *sources, Point *targets, hpx_addr_t rwaddr) {
  RankWise<DualTree> global_tree{rwaddr};
  auto local_tree = global_tree.here();

  // Note: all the pointers are local to the calling rank.
  int first = local_tree->first(rank);
  int last = local_tree->last(rank);
  int range = last - first + 1;
  int send_ns = 0, send_nt = 0;
  for (int i = first; i <= last; ++i) {
    send_ns += count_s[i];
    send_nt += count_t[i];
  }

  // Parcel message length
  size_t bytes = sizeof(hpx_addr_t) + sizeof(int) * 4;
  if (send_ns) {
    bytes += sizeof(int) * range + sizeof(Point) * send_ns;
  }
  if (send_nt) {
    bytes += sizeof(int) * range + sizeof(Point) * send_nt;
  }

  // Acquire parcel
  hpx_parcel_t *p = hpx_parcel_acquire(nullptr, bytes);
  void *data = hpx_parcel_get_data(p);
  hpx_addr_t *rwarg = static_cast<hpx_addr_t *>(data);
  *rwarg = rwaddr;
  int *meta = reinterpret_cast<int *>(
                  static_cast<char *>(data) + sizeof(hpx_addr_t));
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

  char *meta_s = static_cast<char *>(data) + sizeof(hpx_addr_t) +
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
           HPX_POINTER, HPX_POINTER, HPX_ADDR);

// This is the action on the other side that receives the partitioned tree.
int recv_node_handler(char *message_buffer, size_t UNUSED) {
  hpx_addr_t *rwdata = reinterpret_cast<hpx_addr_t *>(message_buffer);
  int *compressed_tree = reinterpret_cast<int *>(
                              message_buffer + sizeof(hpx_addr_t));
  RankWise<DualTree> global_tree{*rwdata};
  auto local_tree = global_tree.here();
  int type = compressed_tree[0];
  int id = compressed_tree[1];
  int n_nodes = compressed_tree[2];
  Node *curr = &(local_tree->unif_grid()[id + type * local_tree->dim3()]);

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
int send_node_handler(Node *n, Point *sorted, int id, int type,
                      hpx_addr_t rwaddr) {
  Node *curr = &n[id];

  RankWise<DualTree> global_tree{rwaddr};
  auto local_tree = global_tree.here();

  // Exclude curr as it is already allocated on remote localities
  int n_nodes = curr->n_descendants() - 1;
  size_t msgsize = sizeof(int) * (3 + n_nodes * 2) + sizeof(hpx_addr_t);
  char *message_buffer = new char[msgsize];
  hpx_addr_t *rwdata = reinterpret_cast<hpx_addr_t *>(message_buffer);
  int *compressed_tree = reinterpret_cast<int *>(
                              message_buffer + sizeof(hpx_addr_t));

  *rwdata = rwaddr;
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
      hpx_parcel_t *p = hpx_parcel_acquire(message_buffer, msgsize);
      hpx_parcel_set_target(p, HPX_THERE(r));
      hpx_parcel_set_action(p, recv_node_action);
      hpx_parcel_send(p, done);
    }
  }

  // Clear local memory once all parcels are sent
  hpx_lco_wait(done);

  delete [] message_buffer;
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, send_node_action, send_node_handler,
           HPX_POINTER, HPX_POINTER, HPX_INT, HPX_INT, HPX_ADDR);


// This appears to be the main action that creates the trees. This will
// organize and call out to the other actions.
// NOTE: One thing is sure, this ought to be factored
int create_dual_tree_handler(hpx_addr_t rwtree, hpx_addr_t sources_gas,
                             hpx_addr_t targets_gas) {
  int rank = hpx_get_my_rank();
  int num_ranks = hpx_get_num_ranks();

  RankWise<DualTree> global_tree{rwtree};
  auto tree = global_tree.here();

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
  int *local_offset_s{nullptr};
  int *local_offset_t{nullptr};
  int *local_count = sort_local_points(&*tree, p_s, n_sources,
                                       p_t, n_targets, &local_offset_s,
                                       &local_offset_t);
  int *local_scount = local_count;
  int *local_tcount = &local_count[tree->dim3()];

  // Compute point distribution
  hpx_lco_get(tree->unif_count(), sizeof(int) * (tree->dim3() * 2),
              tree->unif_count_value());
  tree->set_distribute(distribute_points(num_ranks, tree->unif_count_value(),
                                         tree->dim3()));
  assert(tree->distribute() != nullptr);


  // Exchange points
  Node *ns = tree->unif_grid();
  Node *nt = &ns[tree->dim3()];

  int *global_offset_s = init_point_exchange(rank, &*tree, local_scount,
                                             local_offset_s, p_s, ns, 's');
  int *global_offset_t = init_point_exchange(rank, &*tree, local_tcount,
                                             local_offset_t, p_t, nt, 't');
  hpx_lco_and_set(tree->unif_done(), HPX_NULL);

  // So this one is pretty simple. It sends those points from this rank
  // going to the other rank in a parcel.
  for (int r = 0; r < num_ranks; ++r) {
    if (r != rank) {
      hpx_call(HPX_HERE, send_points_action, HPX_NULL, &r, &local_scount,
               &local_tcount, &local_offset_s, &local_offset_t,
               &p_s, &p_t, &rwtree);
    }
  }

  hpx_addr_t dual_tree_complete = hpx_lco_and_new(2 * tree->dim3());

  for (int r = 0; r < num_ranks; ++r) {
    int first = tree->first(r);
    int last = tree->last(r);
    int s{0}, t{1};

    if (r == rank) {
      for (int i = first; i <= last; ++i) {
        if (tree->unif_count_value()[i] == 0) {
          hpx_lco_and_set(dual_tree_complete, HPX_NULL);
        } else {
          Point *arg = tree->sorted_src();
          hpx_call_when_with_continuation(ns[i].complete(), HPX_HERE,
                                          send_node_action, dual_tree_complete,
                                          hpx_lco_set_action, &ns, &arg,
                                          &i, &s, &rwtree);
        }

        if (tree->unif_count_value()[i + tree->dim3()] == 0) {
          hpx_lco_and_set(dual_tree_complete, HPX_NULL);
        } else {
          Point *arg = tree->sorted_tar();
          hpx_call_when_with_continuation(nt[i].complete(), HPX_HERE,
                                          send_node_action, dual_tree_complete,
                                          hpx_lco_set_action, &nt, &arg,
                                          &i, &t, &rwtree);
        }
      }
    } else {
      for (int i = first; i <= last; ++i) {
        if (tree->unif_count_value()[i] == 0) {
          hpx_lco_and_set(dual_tree_complete, HPX_NULL);
        } else {
          hpx_call_when(ns[i].complete(), dual_tree_complete,
                        hpx_lco_set_action, HPX_NULL, NULL, 0);
        }

        if (tree->unif_count_value()[i + tree->dim3()] == 0) {
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
  delete [] local_offset_s;
  delete [] local_offset_t;
  delete [] global_offset_s;
  delete [] global_offset_t;

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, create_dual_tree_action, create_dual_tree_handler,
           HPX_ADDR, HPX_ADDR, HPX_ADDR);


/////////////////////////////////////////////////////////////////////
// Finalize the partition work - that is, clean up
/////////////////////////////////////////////////////////////////////

int finalize_partition_handler(hpx_addr_t rwtree) {
  RankWise<DualTree> global_tree{rwtree};
  auto tree = global_tree.here();
  tree->clear_data();
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE, finalize_partition_action,
           finalize_partition_handler, HPX_ADDR);


/////////////////////////////////////////////////////////////////////
// DualTree utility stuff
/////////////////////////////////////////////////////////////////////

// This should be called from inside HPX-5.
//
// Also, we want this to return as soon as the work is started. In this way,
// we can do whatever overlap is possible. Then there should be some interface
// to be sure it is complete or something. Possibly even another version of
// create that is create_sync. This means whatever overlap stuff we have going
// will have to be saved in the rankwise data.
RankWise<DualTree> dual_tree_create(int threshold,
                                    Array<Point> sources,
                                    Array<Point> targets) {
  hpx_addr_t domain_geometry = compute_domain_geometry(sources, targets);
  RankWise<DualTree> retval = setup_basic_data(threshold, domain_geometry);
  hpx_lco_delete_sync(domain_geometry);
  return retval;
}


hpx_addr_t dual_tree_partition(RankWise<DualTree> global_tree,
                               Array<Point> sources, Array<Point> targets) {
  hpx_addr_t retval = hpx_lco_future_new(0);
  assert(retval != HPX_NULL);

  hpx_addr_t tree_gas = global_tree.data();
  hpx_addr_t source_gas = sources.data();
  hpx_addr_t target_gas = targets.data();
  hpx_bcast_lsync(create_dual_tree_action, retval,
                  &tree_gas, &source_gas, &target_gas);

  return retval;
}


// This should be called from inside HPX-5
void dual_tree_destroy(RankWise<DualTree> global_tree) {
  hpx_addr_t rwtree = global_tree.data();
  hpx_bcast_rsync(finalize_partition_action, &rwtree);

  auto tree = global_tree.here();
  hpx_lco_delete_sync(tree->unif_count());
}


/////////////////////////////////////////////////////////////////////
// DualTree interface stuff
/////////////////////////////////////////////////////////////////////

void DualTree::clear_data() {
  int rank = hpx_get_my_rank();

  delete [] unif_count_value_;
  hpx_lco_delete_sync(unif_done_);
  delete [] sorted_src_;
  delete [] sorted_tar_;

  int b = first(rank);
  int e = last(rank);

  // NOTE: The difference here is that there are two different allocation
  // schemes for the nodes.

  // TODO: Add some parallelism here. Otherwise, this will take a long time.
  for (int i = 0; i < b; ++i) {
    Node *curr = &unif_grid_[i];
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

  for (int i = b; i <= e; ++i) {
    Node *curr = &unif_grid_[i];
    hpx_lco_delete_sync(curr->sema());
    hpx_lco_delete_sync(curr->complete());

    for (int j = 0; j < 8; ++j) {
      Node *child = curr->child(j);
      if (child) {
        child->destroy(false);
      }
    }
  }

  for (int i = e + 1; i < dim3_; ++i) {
    Node *curr = &unif_grid_[i];
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

  for (int i = 0; i < b; ++i) {
    Node *curr = &unif_grid_[i + dim3_];
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

  for (int i = b; i <= e; ++i) {
    Node *curr = &unif_grid_[i + dim3_];
    hpx_lco_delete_sync(curr->sema());
    hpx_lco_delete_sync(curr->complete());

    for (int j = 0; j < 8; ++j) {
      Node *child = curr->child(j);
      if (child) {
        child->destroy(false);
      }
    }
  }

  for (int i = e + 1; i < dim3_; ++i) {
    Node *curr = &unif_grid_[i + dim3_];
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

  delete [] unif_grid_;
  delete [] distribute_;
}


/////////////////////////////////////////////////////////////////////
// Node interface stuff
/////////////////////////////////////////////////////////////////////

void Node::partition(Point *p, int threshold, DomainGeometry *geo) {
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
    double h = geo->size() / pow(2, idx_.level());
    double center_x = geo->low().x() + (idx_.x() + 0.5) * h;
    double center_y = geo->low().y() + (idx_.y() + 0.5) * h;
    double center_z = geo->low().z() + (idx_.z() + 0.5) * h;

    // Sort the particles among the children
    Point *splits[9]{};
    splits[0] = &p[first_];
    splits[8] = &p[last_ + 1];

    auto z_comp = [&center_z](Point &a) {
      return a.z() < center_z;
    };
    splits[4] = std::partition(splits[0], splits[8], z_comp);

    auto y_comp = [&center_y](Point &a) {
      return a.y() < center_y;
    };
    splits[2] = std::partition(splits[0], splits[4], y_comp);
    splits[6] = std::partition(splits[4], splits[8], y_comp);

    auto x_comp = [&center_x](Point &a) {
      return a.x() < center_x;
    };
    splits[1] = std::partition(splits[0], splits[2], x_comp);
    splits[3] = std::partition(splits[2], splits[4], x_comp);
    splits[5] = std::partition(splits[4], splits[6], x_comp);
    splits[7] = std::partition(splits[6], splits[8], x_comp);

    // Perform some counting
    int stat[8]{};
    int offset[8]{};

    offset[0] = first_;
    stat[0] = splits[1] - splits[0];
    for (int i = 1; i < 8; ++i) {
      stat[i] = splits[i + 1] - splits[i];
      offset[i] = offset[i - 1] + stat[i - 1];
    }

    // Create child nodes
    for (int i = 0; i < 8; ++i) {
      if (stat[i]) {
        Node *child = new Node{idx_.child(i), offset[i],
                               offset[i] + stat[i] - 1, this};
        child_[i] = child;

        hpx_call(HPX_HERE, partition_node_action, HPX_NULL,
                 &child, &p, &geo, &threshold);
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
  if (sema_ != HPX_NULL) {
    hpx_lco_delete_sync(sema_);
  }
  if (complete_ != HPX_NULL) {
    hpx_lco_delete_sync(complete_);
  }
  if (allocated_in_array == false) {
    delete this;
  }
}

