#include <cstdlib>
#include <cmath>
#include <cstdint> 
#include <cstring>
#include <iostream>
#include "tree.h"

uint64_t split(unsigned k) {  
  uint64_t split = k & 0x1fffff; 
  split = (split | split << 32) & 0x1f00000000ffff;
  split = (split | split << 16) & 0x1f0000ff0000ff;
  split = (split | split << 8)  & 0x100f00f00f00f00f;
  split = (split | split << 4)  & 0x10c30c30c30c30c3;
  split = (split | split << 2)  & 0x1249249249249249;
  return split;
}

uint64_t morton_key(const index_t *index) {
  uint64_t key = 0;
  key |= split(x) | split(y) << 1 | split(z) << 2;
  return key;
}

int allocate_points_handler(int nsrc, int ntar, char datatype, 
                            hpx_addr_t *src_handle, 
                            hpx_addr_t *tar_handle) {
  int num_ranks = hpx_get_num_ranks(); 
  hpx_addr_t src = hpx_gas_calloc_cyclic(num_ranks, sizeof(Point) * nsrc, 0); 
  hpx_addr_t tar = hpx_gas_calloc_cyclic(num_ranks, sizeof(Point) * ntar, 0); 
  *src_handle = src; 
  *tar_handle = tar; 
  hpx_bcast_rsync(set_points_action, &src, &tar, &nsrc, &ntar, &datatype); 
  hpx_exit(HPX_SUCCESS); 
}
HPX_ACTION(HPX_DEFAULT, 0, allocate_points_action, allocate_points_handler, 
           HPX_INT, HPX_INT, HPX_CHAR, HPX_POINTER, HPX_POINTER); 


int set_points_handler(hpx_addr_t src, hpx_addr_t tar, int nsrc, int ntar, 
                       char datatype) {
  int rank = hpx_get_my_rank(); 
  hpx_addr_t curr_src = hpx_addr_add(src, sizeof(Point) * nsrc * rank, 
                                     sizeof(Point) * nsrc); 
  hpx_addr_t curr_tar = hpx_addr_add(tar, sizeof(Point) * ntar * rank,
                                     sizeof(Point) * ntar); 
  srand(rank); 
  Point *s{nullptr}, *t{nullptr}; 
  
  if (!hpx_gas_try_pin(curr_src, (void **)&s) || 
      !hpx_gas_try_pin(curr_tar, (void **)&t)) 
    return HPX_ERROR; 

  if (datatype == 'c') {
    for (int i = 0; i < nsrc; ++i) {
      double x = 1.0 * rand() / RAND_MAX - 0.5; 
      double y = 1.0 * rand() / RAND_MAX - 0.5; 
      double z = 1.0 * rand() / RAND_MAX - 0.5; 
      s[i] = Point{x, y, z}; 
    }

    for (int i = 0; i < ntargets; ++i) {
      double x = 1.0 * rand() / RAND_MAX - 0.5; 
      double y = 1.0 * rand() / RAND_MAX - 0.5; 
      double z = 1.0 * rand() / RAND_MAX - 0.5; 
      t[i] = Point{x, y, z}; 
    }
  } else {
    for (int i = 0; i < nsources; ++i) {
      double theta = 1.0 * rand() / RAND_MAX * M_PI_2; 
      double phi = 1.0 * rand() / RAND_MAX * M_PI * 2; 
      double x = sin(theta) * cos(phi); 
      double y = sin(theta) * sin(phi); 
      double z = cos(theta); 
      s[i] = Point{x, y, z}; 
    }

    for (int i = 0; i < ntargets; ++i) {
      double theta = 1.0 * rand() / RAND_MAX * M_PI_2; 
      double phi = 1.0 * rand() / RAND_MAX * M_PI * 2; 
      double x = sin(theta) * cos(phi); 
      double y = sin(theta) * sin(phi); 
      double z = cos(theta); 
      t[i] = Point{x, y, z}; 
    }
  }

  hpx_gas_unpin(curr_src); 
  hpx_gas_unpin(curr_tar); 
  return HPX_SUCCESS; 
}
HPX_ACTION(HPX_DEFAULT, 0, set_points_action, set_points_handler, 
           HPX_ADDR, HPX_ADDR, HPX_INT, HPX_INT, HPX_CHAR);

int delete_points_handler(hpx_addr_t src, hpx_addr_t tar) { 
  hpx_gas_free_sync(src); 
  hpx_gas_free_sync(tar); 
  hpx_exit(HPX_SUCCESS); 
} 
HPX_ACTION(HPX_DEFAULT, 0, delete_points_action, delete_points_handler, 
           HPX_ADDR, HPX_ADDR); 


int partition_points_handler(int nsrc, hpx_addr_t src, 
                             int ntar, hpx_addr_t tar, int threshold) {
  int num_ranks = hpx_get_num_ranks(); 

  // Determine domain geometry
  hpx_addr_t domain_geometry = 
    hpx_lco_user_new(sizeof(int) + sizeof(double) * 6, 
                     domain_geometry_init_action, 
                     domain_geometry_op_action, 
                     domain_geometry_predicate_action, 
                     &num_ranks, sizeof(num_ranks)); 

  hpx_bcast_lsync(set_domain_geometry_action, HPX_NULL, 
                  &nsrc, &src, &ntar, &tar, &domain_geometry); 

  void *temp{nullptr}; 
  hpx_lco_getref{domain_geometr, 1, &temp); 
  double *var = reinterpret_cast<double *>(static_cast<char *>(temp) + 
                                           sizeof(int)); 
  double size_x = var[1] - var[0]; 
  double size_y = var[3] - var[2]; 
  double size_z = var[5] - var[4]; 
  double size = fmax(size_x, fmax(size_y, size_z)); 
  double corner_x = (var[1] + var[0] - size) / 2; 
  double corner_y = (var[3] + var[2] - size) / 2; 
  double corner_z = (var[5] + var[4] - size) / 2; 

  // Choose uniform partition level such that the number of grids is no less
  // than the number of ranks 
  unif_level = ceil(log(num_ranks) / log(8)); 
  int n_grids = pow(8, unif_level); 

  // Cyclic allocation on rank 0
  unif_count = hpx_gas_calloc_cyclic(num_ranks, sizeof(hpx_addr_t), 0); 
  unif_grid = hpx_gas_calloc_cyclic(num_ranks, sizeof(Node) * n_grids * 2, 0);
  unif_done = hpx_gas_calloc_cyclic(num_ranks, sizeof(hpx_addr_t), 0); 
  sorted_src = hpx_gas_calloc_cyclic(num_ranks, sizeof(Point *), 0); 
  sorted_tar = hpx_gas_calloc_cyclic(num_ranks, sizeof(Point *), 0); 

  hpx_bcast_rsync(init_partition_action, &unif_count, &unif_grid, &unif_done, 
                  &sorted_src, &sorted_tar, &unif_level, &threshold, 
                  &corner_x, &corner_y, &corner_z, &size); 

  hpx_bcast_rsync(create_dual_tree_action, &src, &nsrc, &tar, &ntar); 

  // Cleanup 
  // unif_count, unif_grid, unif_done, sorted_src, sorted_tar

  hpx_exit(HPX_SUCCESS); 
}
HPX_ACTION(HPX_DEFAULT, 0, partition_points_action, partition_points_handler, 
           HPX_INT, HPX_ADDR, HPX_INT, HPX_ADDR, HPX_INT); 

int init_partition_handler(hpx_addr_t count, hpx_addr_t grid, hpx_addr_t done,
                           hpx_addr_t src, hpx_addr_t tar, int level, 
                           int limit, double cx, double cy, double cz, 
                           double sz) {
  int rank = hpx_get_my_rank(); 
  int num_ranks = hpx_get_num_ranks(); 
  int dim = pow(2, level), dim3 = pow(8, level); 

  if (rank) {
    unif_count = count; 
    unif_grid = grid; 
    unif_done = done; 
    sorted_src = src; 
    sorted_tar = tar; 
    unif_level = level;
    threshold = limit; 
    corner_x = cx; 
    corner_y = cy; 
    corner_z = cz; 
    size = sz; 
  }

  // Setup unif_count LCO 
  hpx_addr_t curr_unif_count = hpx_addr_add(unif_count, 
                                            sizeof(hpx_addr_t) * rank, 
                                            sizeof(hpx_addr_t)); 
  hpx_addr_t *user_lco{nullptr}; 
  if (!hpx_gas_try_pin(curr_unif_count, (void **)&user_lco))
    return HPX_ERROR; 

  int init[2] = {num_ranks, dim3 * 2}; 
  *user_lco = hpx_lco_user_new(sizeof(int) * (2 + dim3 * 2), 
                               unif_grid_count_init_action, 
                               unif_grid_count_op_action, 
                               unif_grid_count_predicate_action,
                               &init, sizeof(init)); 
  hpx_gas_unpin(curr_unif_count); 
 
  // Setup unif_grid
  hpx_addr_t curr_unif_grid = hpx_addr_add(unif_grid, 
                                           sizeof(Node) * dim3 * 2 * rank, 
                                           sizeof(Node) * dim3 * 2); 
  
  Node *n{nullptr}, *s{nullptr}, *t{nullptr}; 
  if (!hpx_gas_try_pin(curr_unif_grid, (void **)&n))
    return HPX_ERROR; 

  s = &n[0]; 
  t = &n[dim3]; 

  int iterator = 0; 
  for (int iz = 0; iz < num_per_dim; ++iz) {
    for (int iy = 0; iy < num_per_dim; ++iy) {
      for (int ix = 0; ix < num_per_dim; ++ix) {
        // The purpose here is to setup Index and LCOs of the Node 
        s[iterator] = Node{Index{level, ix, iy, iz}}; 
        t[iterator] = Node{Index{level, ix, iy, iz}}; 
        iterator++; 
      }
    }
  }

  hpx_gas_unpin(curr_unif_grid); 

  // Setup unif_done
  hpx_addr_t curr_unif_done = hpx_addr_add(unif_done, 
                                           sizeof(hpx_addr_t) * rank, 
                                           sizeof(hpx_addr_t)); 
  hpx_addr_t *gate{nullptr}; 
  if (!hpx_gas_try_pin(curr_unif_done, (void **)&gate))
    return HPX_ERROR; 

  *gate = hpx_lco_and_new(1); 
  hpx_gas_unpin(curr_unif_done); 

  return HPX_SUCCESS; 
} 
HPX_ACTION(DEFAULT, 0, init_partition_action, init_partition_handler, 
           HPX_ADDR, HPX_ADDR, HPX_ADDR, HPX_ADDR, HPX_ADDR, HPX_INT, HPX_INT, 
           HPX_DOUBLE, HPX_DOUBLE, HPX_DOUBLE, HPX_DOUBLE); 


int create_dual_tree_handler(hpx_addr_t src, int nsrc, 
                             hpx_addr_t tar, int ntar) {
  int num_ranks = hpx_get_num_ranks(); 
  int rank = hpx_get_my_rank(); 

  hpx_addr_t curr_src = hpx_addr_add(src, sizeof(Point) * nsrc * rank, 
                                     sizeof(Point) * nsrc); 
  hpx_addr_t curr_tar = hpx_addr_add(tar, sizeof(Point) * ntar * rank, 
                                     sizeof(Point) * ntar); 
  Point *s{nullptr}, *t{nullptr}; 
  if (!hpx_gas_try_pin(curr_src, (void **)&s) || 
      !hpx_gas_try_pin(curr_tar, (void **)&t)) 
    return HPX_ERROR; 

  // Assign points to grid 
  int dim = pow(2, unif_level),  dim3 = pow(8, unif_level); 
  int *gid_of_src = new int[nsrc]; 
  int *gid_of_tar = new int[ntar]; 
  int *local_cnt = new int[dim3 * 2]; 
  int *local_scnt = &local_cnt[0]; 
  int *local_tcnt = &local_cnt[dim3]; 
  double scale = 1.0 / size; 

  for (int i = 0; i < nsrc; ++i) {
    Point *j = &s[i]; 
    int xid = dim * (j->x() - corner_x) * scale; 
    int yid = dim * (j->y() - corner_y) * scale; 
    int zid = dim * (j->z() - corner_z) * scale; 
    int gid = morton_key(xid, yid, zid); 
    gid_of_src[i] = gid; 
    local_scnt[gid]++; 
  }

  for (int i = 0; i < ntar; ++i) {
    Point *j = &t[i]; 
    int xid = dim * (j->x() - corner_x) * scale; 
    int yid = dim * (j->y() - corner_y) * scale; 
    int zid = dim * (j->z() - corner_z) * scale; 
    int gid = morton_key(xid, yid, zid); 
    gid_of_tar[i] = gid; 
    local_tcnt[gid]++;
  }

  // Exchange counts
  for (int r = 0; i < num_ranks; ++r) {
    hpx_parcel_t *p = hpx_parcel_acquire(local_cnt, sizeof(int) * dim3 * 2); 
    hpx_parcel_set_target(p, HPX_THERE(r)); 
    hpx_parcel_set_action(p, exchange_count_action); 
    hpx_parcel_send(p, HPX_NULL); 
  }

  // Put points in the same grid together  
  int *local_offset_s = new int[dim3]; 
  int *local_offset_t = new int[dim3]; 
  Point *temp_s = new Point[nsrc]; 
  Point *temp_t = new Point[ntar]; 
  int *assigned = new int[dim3]; 

  for (int i = 1; i < dim3; ++i) {
    local_offset_s[i] = local_offset_s[i - 1] + local_scnt[i - 1]; 
    local_offset_t[i] = local_offset_t[i - 1] + local_tcnt[i - 1]; 
  }

  for (int i = 0; i < nsrc; ++i) {
    int gid = gid_of_src[i]; 
    int cursor = local_offset_s[gid] + assigned[gid]; 
    temp_s[cursor] = s[i]; 
    assigned[gid]++;
  }

  memset(assigned, 0, dim3); 

  for (int i = 0; i < ntar; ++i) {
    int gid = gid_of_tar[i]; 
    int cursor = local_offset_t[gid] + assigned[gid]; 
    temp_t[cursor] = t[i]; 
    assigned[gid]++;
  }

  delete [] gid_of_src; 
  delete [] gid_of_tar; 
  delete [] assigned; 
    
  hpx_gas_unpin(curr_src); 
  hpx_gas_unpin(curr_tar); 

  // Compute point distribution
  hpx_addr_t curr_unif_count = hpx_addr_add(unif_count, 
                                            sizeof(hpx_addr_t) * rank, 
                                            sizeof(hpx_addr_t)); 
  hpx_addr_t *user_lco{nullptr}; 
  if (!hpx_gas_try_pin(curr_unif_count, (void **)&user_lco))
    return HPX_ERROR; 
  
  void *user_lco_buffer{nullptr}; 
  int *global_cnt{nullptr}; 
  hpx_lco_getref(*user_lco, 1, &user_lco_buffer); 
  global_cnt = reinterpret_cast<int *>(static_cast<char *>(user_lco_buffer) 
                                       + sizeof(int) * 2); 
  int *dist = distribute_points(num_ranks, global_cnt, dim3); 

  // Exchange points
  hpx_addr_t curr_sorted_src = hpx_addr_add(sorted_src, 
                                            sizeof(Point *) * rank, 
                                            sizeof(Point *)); 
  hpx_addr_t curr_sorted_tar = hpx_addr_add(sorted_tar, 
                                            sizeof(Point *) * rank, 
                                            sizeof(Point *)); 
  hpx_addr_t curr_unif_done = hpx_addr_add(unif_done, 
                                           sizeof(hpx_addr_t) * rank, 
                                           sizeof(hpx_addr_t)); 
  hpx_addr_t curr_unif_grid = hpx_addr_add(unif_grid, 
                                           sizeof(Node) * dim3 * 2 * rank, 
                                           sizeof(Node) * dim3 * 2); 

  Point **ss{nullptr}, **st{nullptr}; 
  Node *n{nullptr}, *ns{nullptr}, *nt{nullptr}; 
  hpx_addr_t *gate{nullptr}; 

  if (!hpx_gas_try_pin(curr_sorted_src, (void **)&ss) ||
      !hpx_gas_try_pin(curr_sorted_tar, (void **)&st) ||
      !hpx_gas_try_pin(curr_unif_done, (void **)&gate) ||
      !hpx_gas_try_pin(curr_unif_grid, (void **)&n))
    return HPX_ERROR; 

  ns = &n[0]; 
  nt = &n[dim3]; 

  // Range of grids owned by the current rank
  int start = (rank == 0 ? 0 : dist[rank - 1] + 1); 
  int end = dist[rank]; 
  int range = end - start + 1; 
  int nss = global_cnt[start], nst = global_cnt[start + dim3]; 
  int *global_offset_s = new int[range]; 
  int *global_offset_t = new int[range]; 

  for (int i = 1; i < range; ++i) {
    nss += global_cnt[i + start]; 
    nst += global_cnt[i + start + dim3]; 
    global_offset_s[i] = global_offset_s[i - 1] + global_cnt[i + start];
    global_offset_t[i] = global_offset_t[i - 1] + global_cnt[i + start]; 
  }

  // Allocate memory to hold points of the assigned grids
  *ss = new Point[nss]; 
  *st = new Point[nst]; 

  // And memory for adaptive partitioning 
  swap_src = new int[nss]; 
  bin_src = new int[nss]; 
  map_src = new int[nss]; 
  swap_tar = new int[nst]; 
  bin_tar = new int[nst]; 
  map_tar = new int[nst]; 

  for (int i = 0; i < nss; ++i) 
    map_src[i] = i; 
  for (int i = 0; i < nst; ++i) 
    map_tar[i] = i; 

  for (int i = start; i <= end; ++i) {
    // Copy local source points
    memcpy(*ss + global_offset_s[i - start], &temp_s[local_offset_s[i]],
           sizeof(Point) * local_scnt[i]); 

    // Here, first_ will be used as an iterator to trace the position to 
    // insert point. When all the points have been inserted, the correct 
    // value of first_ can be recovered from last_ and global_cnt. 
    ns[i].set_first(global_offset_s[i - start] + local_scnt[i]); 
    
    if (i < end) {
      ns[i].set_last(global_offset_s[i + 1 - start]); 
    } else {
      ns[i].set_last(nss - 1); 
    }

    // Copy local target points 
    memcpy(*st + global_offset_t[i - start], &temp_t[local_offset_t[i]], 
           sizeof(Point) * local_tcnt[i]); 

    nt[i].set_first(global_offset_t[i - start] + local_tcnt[i]); 
    
    if (i < end) {
      nt[i].set_last(global_offset_t[i + 1 - start]); 
    } else {
      nt[i].set_last(nst - 1);
    }
  }

  // The current rank is ready to receive 
  hpx_lco_and_set(*gate, HPX_NULL); 

  // Send points to the other ranks 
  for (int r = 0; r < num_ranks; ++r) {
    if (r != rank) {
      // Range of grids assigned to rank r
      int start = (r == 0 ? 0 : dist[r - 1] + 1); 
      int end = dist[r]; 
      int range = end - start + 1; 

      int send_ns = local_offset_s[end] - local_offset_s[start] + 1; 
      int send_nt = local_offset_t[end] - local_offset_t[start] + 1; 

      // Parcel message size 
      size_t psize = sizeof(int) * (4 + 2 * range) + 
        sizeof(Point) * (send_ns + send_nt); 
      hpx_parcel_t *p = hpx_parcel_acquire(NULL, psize); 

      void *data = hpx_parcel_get_data(p); 

      int *meta = reinterpret_cast<int *>(static_cast<char *>(data)); 
      meta[0] = start;  
      meta[1] = end; 
      meta[2] = send_ns; 
      meta[3] = send_nt; 
      int *meta_offset_s = &meta[4]; 
      int *meta_offset_t = &meta[4 + range]; 

      for (int i = start; i <= end; ++i) {
        meta_offset_s[i - start] = local_offset_s[i] - local_offset_s[start]; 
        meta_offset_t[i - start] = local_offset_t[i] - local_offset_t[start]; 
      }

      char *meta_s = static_cast<char *>(data) + sizeof(int) * (4 + range * 2); 
      char *meta_t = meta_s + sizeof(Point) * send_ns; 

      memcpy(meta_s, &temp_s[start], send_ns); 
      memcpy(meta_t, &temp_t[start], send_nt); 

      hpx_parcel_set_target(p, HPX_THERE(r)); 
      hpx_parcel_set_action(p, exchange_points_action); 
      hpx_parcel_send(p, HPX_NULL);
    }
  }

  // Exchange source and target trees
  hpx_addr_t dual_tree_complete = hpx_lco_and_new(2 * dim3); 

  for (int i = 0; i < start; ++i) {
    hpx_call_when(ns[i].complete(), dual_tree_complete, hpx_lco_set_action, 
                  HPX_NULL, NULL, 0); 
    hpx_call_when(nt[i].complete(), dual_tree_complete, hpx_lco_set_action, 
                  HPX_NULL, NULL, 0); 
  }

  for (int i = end + 1; i < dim3; ++i) {
    hpx_call_when(ns[i].complete(), dual_tree_complete, hpx_lco_set_action, 
                  HPX_NULL, NULL, 0); 
    hpx_call_when(nt[i].complete(), dual_tree_complete, hpx_lco_set_action, 
                  HPX_NULL, NULL, 0); 
  }    

  for (int i = start; i <= end; ++i) {
    int s{0}, t{1}; 
    hpx_call_when_with_continuation(ns[i].complete(), HPX_HERE, 
                                    send_node_action, 
                                    dual_tree_complete, hpx_lco_set_action, 
                                    &ns, ss, &i, &s); 

    hpx_call_when_with_continuation(nt[i].complete(), HPX_HERE, 
                                    send_node_action, 
                                    dual_tree_complete, hpx_lco_set_action, 
                                    &nt, st, &i, &t); 
  }
      
  hpx_lco_wait(dual_tree_complete); 
  hpx_lco_delete_sync(dual_tree_complete); 

  hpx_lco_release(*user_lco, user_lco_buffer); 
  hpx_gas_unpin(curr_unif_count); 
  hpx_gas_unpin(curr_sorted_src); 
  hpx_gas_unpin(curr_sorted_tar); 
  hpx_gas_unpin(curr_unif_done); 
  hpx_gas_unpin(curr_unif_grid); 

  delete [] local_cnt; 
  delete [] local_offset_s; 
  delete [] local_offset_t; 
  delete [] temp_s; 
  delete [] temp_t; 
  delete [] swap_src; 
  delete [] bin_src; 
  delete [] map_src; 
  delete [] swap_tar; 
  delete [] bin_tar; 
  delete [] map_tar; 

  return HPX_SUCCESS; 
}
HPX_ACTION(HPX_DEFAULT, 0, create_dual_tree_action, create_dual_tree_handler, 
           HPX_ADDR, HPX_INT, HPX_ADDR, HPX_INT, HPX_DOUBLE, HPX_DOUBLE, 
           HPX_DOUBLE, HPX_DOUBLE); 

int send_node_handler(Node *n, Point *p, int id, int type) {
  int first = n[id].first(); 
  int last = n[id].last(); 
  int range = last - first = 1; 

  // Rearrange points
  int *map = (type ? map_tar : map_src); 
  Point *temp = new Point[range]; 
  for (int i = first; i <= last; ++i) 
    temp[i - first] = p[map[i]]; 
  memcpy(p + first, temp, sizeof(Point) * range); 
  delete [] temp; 

  // Exclude n as it is already allocated on remote localities
  int n_nodes = n->descendants() - 1; 
  int *compressed_tree = new int[3 + n_nodes * 2]; 

  compressed_tree[0] = type; // target is 0, source is 1
  compressed_tree[1] = id; // where to merge the tree 
  compressed_tree[2] = n_nodes; // # of nodes 
  int *branch = &compressed_tree[3]; 
  int *tree = &compressed_tree[3 + n_nodes]; 
  int curr = 0; 
  n->compress(branch, tree, -1, curr); 
  
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
           HPX_POINTER_T, HPX_POINTER_T, HPX_INT, HPX_INT); 

int recv_node(void *args, size_t size) {
  int *compressed_tree = reinterpret_cast<int *>(static_cast<char *>(args)); 

  int rank = hpx_get_my_rank(); 
  int dim3 = pow(8, unif_level); 
  hpx_addr_t curr_unif_grid = hpx_addr_add(unif_grid, 
                                           sizeof(Node) * dim3 * 2 * rank, 
                                           sizeof(Node) * dim3 * 2); 
  Node *n{nullptr}, *curr{nullptr}; 
  if (!hpx_gas_try_pin(curr_unif_grid, (void **)&n))
    return HPX_ERROR; 
  
  int type = compressed_tree[0]; 
  int id = compressed_tree[1]; 
  int n_nodes = compressed_tree[2]; 
  const int *branch = &compressed_tree[3]; 
  const int *tree = &compressed_tree[3 + n_nodes]; 

  curr = &n[id + type * n_nodes]; 
  curr->extract(branch, tree, n_nodes); 
  
  hpx_lco_and_set_num(curr->complete(), 8, HPX_NULL);   
  hpx_gas_unpin(curr_unif_grid);

  return HPX_SUCCESS; 
} 
HPX_ACTION(HPX_DEFAULT, 0, recv_node_action, recv_node_handler, 
           HPX_POINTER_T, HPX_SIZE_T); 

int exchange_point_handler(void *args, size_t size) {
  int rank = hpx_get_my_rank(); 
  hpx_addr_t curr_unif_done = hpx_addr_t(unif_done, 
                                         sizeof(hpx_addr_t) * rank, 
                                         sizeof(hpx_addr_t)); 
  hpx_addr_t *gate{nullptr}; 
  if (!hpx_gas_try_pin(curr_unif_done, (void **)&gate)) 
    return HPX_ERROR; 

  // Wait until the buffer for exchanging points has been allocated 
  hpx_lco_wait(*gate); 
  hpx_gas_unpin(curr_unif_done); 

  // Now proceed to merge the incoming points 
  int dim = pow(2, unif_level), dim3 = pow(8, unif_level); 

  hpx_addr_t curr_unif_count = hpx_addr_t(unif_count, 
                                          sizeof(hpx_addr_t) * rank, 
                                          sizeof(hpx_addr_t)); 
  hpx_addr_t curr_sorted_src = hpx_addr_t(sorted_src, 
                                          sizeof(Point *) * rank, 
                                          sizeof(Point *)); 
  hpx_addr_t curr_sorted_tar = hpx_addr_t(sorted_tar, 
                                          sizeof(Point *) * rank, 
                                          sizeof(Point *)); 
  hpx_addr_t curr_unif_grid = hpx_addr_t(unif_grid, 
                                         sizeof(Node) * dim3 * 2 * rank, 
                                         sizeof(Node) * dim3 * 2); 
  Point **s{nullptr}, **t{nullptr}; 
  Node *n{nullptr}, *ns{nullptr}, *nt{nullptr}; 
  hpx_addr_t *user_lco{nullptr}; 
  if (!hpx_gas_try_pin(curr_unif_count, (void **)&user_lco) ||
      !hpx_gas_try_pin(curr_sorted_src, (void **)&s) ||
      !hpx_gas_try_pin(curr_sorted_tar, (void **)&t) ||
      !hpx_gas_try_pin(curr_unif_grid, (void **)&n))
    return HPX_ERROR; 

  void *user_lco_buffer{nullptr}; 
  hpx_lco_getref(*user_lco, 1, &user_lco_buffer); 
  int *count = reinterpret_cast<int *>(static_cast<char *>(user_lco_buffer) 
                                       + sizeof(int) * 2); 
  int *count_s = count; 
  int *count_t = &count[dim3]; 

  // Process incoming message 
  int *meta = reinterpret_cast<int *>(static_cast<char *>(args)); 
  int start = meta[0];
  int end = meta[1];
  int range = end - start + 1; 
  int recv_ns = meta[2]; 
  int recv_nt = meta[3]; 
  int *offset_s = &meta[4]; 
  int *offset_t = &meta[4 + range]; 
  
  Point *recv_s = reinterpret_cast<Point *>(static_cast<char *>(args) + 
                                            sizeof(int) * (4 + range * 2)); 
  Point *recv_t = recv_s + recv_ns; 

  for (int i = start; i <= end; ++i) {
    Node *curr_s = &ns[i]; 
    Node *curr_t = &nt[i]; 

    int incoming_first_s = offset_s[i - start]; 
    int incoming_last_s = (i == end ? recv_ns - 1 : 
                           offset_s[i + 1 - start] - 1); 
    int incoming_ns = incoming_last_s - incoming_first_s + 1; 

    hpx_lco_sema_p(curr_s->sema()); 
    int curr_s_first = curr_s->first(); 
    memcpy(*s + curr_s_first, &recv_s[incoming_first_s], 
           sizeof(Point) * incoming_ns); 
    curr_s->set_first(curr_s_first + incoming_ns); 

    if (curr_s_first + incoming_ns > curr_s->last()) {
      // This grid is ready for adaptive partition
      int first = curr_s->last() + 1 - count_s[i]; 
      curr_s->set_first(first); 
      
      hpx_call(HPX_HERE, partition_node_action, HPX_NULL, 
               &curr_s, s, &swap_src, &bin_src, &map_src); 
    }
    hpx_lco_sema_v(curr_s->sema()); 
  
    int incoming_first_t = offset_t[i - start]; 
    int incoming_last_t = (i == end ? recv_nt - 1:
                           offset_t[i + 1 - start] - 1); 
    int incoming_nt = incoming_last_t - incoming_first_t + 1; 

    hpx_lco_sema_p(curr_t->sema()); 
    int curr_t_first = curr_t->first(); 
    memcpy(*t + curr_t_first, &recv_t[incoming_first_t], 
           sizeof(Point) * incoming_nt); 
    curr_t->set_first(curr_t_first + incoming_nt); 

    if (curr_t_first + incoming_nt > curr_t>last()) {
      // This grid is ready for adaptive partition 
      int first = curr_t->last() + 1 - count_t[i]; 
      curr_t->set_first(first); 

      hpx_call(HPX_HERE, partition_node_action, HPX_NULL, 
               &curr_t, t, &swap_tar, &bin_tar, &map_tar); 
    }
    hpx_lco_sema_v(curr_t->sema()); 
  }

  hpx_lco_release(*user_lco, user_lco_buffer); 
  hpx_gas_unpin(curr_unif_count); 
  hpx_gas_unpin(curr_sorted_src); 
  hpx_gas_unpin(curr_sorted_tar); 
  hpx_gas_unpin(curr_unif_grid); 

  return HPX_SUCCESS; 
} 
HPX_ACTION(HPX_DEFAULT, 0, exchange_point_action, exchange_point_handler, 
           HPX_POINTER, HPX_SIZE_T); 

int partition_node_handler(Node *n, Point *p, int *swap, int *bin, int *map) {
  n->partition(p, swap, bin, map, threshold, corner_x, corner_y, 
               corner_z, size); 
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, partition_node_action, partition_node_handler, 
           HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_POINTER); 

void unif_grid_count_init_handler(void *data, const size_t size, 
                                  int *init, size_t init_size) {

  int *input = reinterpret_cast<int *>(static_cast<char *>(data)); 
  
  input[0] = init[0]; // # of inputs; 
  input[1] = init[1]; // # of terms 

  for (int i = 0; i < init[1]; ++i) 
    input[i + 2] = 0; 
} 
HPX_ACTION(HPX_FUNCTION, 0, unif_grid_count_init_action, 
           unif_grid_count_init_handler, HPX_POINTER, HPX_SIZE_T, 
           HPX_POINTER, HPX_SIZE_T); 

void unif_grid_count_op_handler(void *data, int *rhs, size_t size) {
  int *input = reinterpret_cast<int *>(static_cast<char *>(data)); 
  int *lhs = &input[2]; 

  input[0]--; 
  int nterms = input[1]; 
  for (int i = 0; i < nterms; ++i) 
    lhs[i] += rhs[i]; 
}
HPX_ACTION(HPX_FUNCTION, 0, unif_grid_count_op_action, 
           unif_grid_count_op_handler, HPX_POINTER, HPX_POINTER, HPX_SIZE_T); 

bool unif_grid_count_predicate_handler(void *data, size_t size) {
  int *input = reinterpret_cast<int *>(static_cast<char *>(data)); 
  return (input[0] == 0); 
}
HPX_ACTION(HPX_FUNCTION, 0, unif_grid_count_predicate_action, 
           unif_grid_count_predicate_handler, HPX_POINTER, HPX_SIZE_T); 

int exchange_count_handler(void *args, size_t size) {
  int rank = hpx_get_my_rank(); 
  hpx_addr_t curr_unif_count = hpx_addr_add(unif_count, 
                                            sizeof(hpx_addr_t) * rank, 
                                            sizeof(hpx_addr_t)); 
  hpx_addr_t *user_lco{nullptr}; 
  if (!hpx_gas_try_pin(curr_unif_count, (void **)&user_lco))
    return HPX_ERROR; 
  
  hpx_lco_set(*user_lco, size, args, HPX_NULL, HPX_NULL); 
  hpx_gas_unpin(curr_unif_count); 
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, exchange_count_action, exchange_count_handler, 
           HPX_POINTER, HPX_SIZE_T); 


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
 
int set_domain_geometry_handler(int nsrc, hpx_addr_t src, 
                                int ntar, hpx_addr_t tar, 
                                hpx_addr_t domain_geometry) {
  int rank = hpx_get_my_rank(); 
  hpx_addr_t curr_src = hpx_addr_add(src, sizeof(Point) * nsrc * rank, 
                                     sizeof(Point) * nsrc); 
  hpx_addr_t curr_tar = hpx_addr_add(tar, sizeof(Point) * ntar * rank, 
                                     sizeof(Point) * ntar); 
  Point *s{nullptr}, *t{nullptr}; 
  if (!hpx_gas_try_pin(curr_src, (void **)&s) || 
      !hpx_gas_try_pin(curr_tar, (void **)&t))
    return HPX_ERROR; 

  double var[6] = {1e50, -1e50, 1e50, -1e50, 1e50, -1e50}; 

  for (int i = 0; i < nsrc; ++i) {
    var[0] = fmin(var[0], s[i].x()); 
    var[1] = fmax(var[1], s[i].x()); 
    var[2] = fmin(var[2], s[i].y()); 
    var[3] = fmax(var[3], s[i].y()); 
    var[4] = fmin(var[4], s[i].z()); 
    var[5] = fmax(var[5], s[i].z()); 
  }

  for (int i = 0; i < ntar; ++i) {
    var[0] = fmin(var[0], tar[i].x()); 
    var[1] = fmax(var[1], tar[i].x()); 
    var[2] = fmin(var[2], tar[i].y()); 
    var[3] = fmax(var[3], tar[i].y()); 
    var[4] = fmin(var[4], tar[i].z()); 
    var[5] = fmax(var[5], tar[i].z()); 
  }
  
  hpx_gas_unpin(curr_src); 
  hpx_gas_unpin(curr_tar); 
  
  hpx_lco_set_lsync(domain_geometry, sizeof(double) * 6, var, HPX_NULL);

  return HPX_SUCCESS; 
}
HPX_ACTION(HPX_DEFAULT, 0, set_domain_geometry_action, 
           set_domain_geometry_handler, HPX_INT, HPX_ADDR, 
           HPX_INT, HPX_ADDR, HPX_ADDR); 


int *distribute_points(int num_ranks, const int *global, int len) {
  int *ret = new int[num_ranks]; 

  int *s = global; // Source counts
  int *t = &global[len]; // Target counts

  int *scan = new int[len]; 
  scan[0] = s[0] + t[0]; 
  for (int i = 1; i < len; ++i) 
    scan[i] = scan[i - 1] + s[i] + t[i]; 

  int q = scan[len - 1] / num_ranks; 
  int r = scan[len - 1] % num_ranks; 

  int rank = 0; 
  int iterator = -1; 
  int bound = q + (r != 0); 

  while (rank < num_ranks) {
    if (rank) 
      bound = scan[iterator] + q + (rank <= r); 
    iterator++; 
    for (int i = iterator; i < len; ++i) {
      if (scan[i] <= bound && scan[i + 1] > bound)
        break;
    }
    ret[rank++] = i; 
    iterator = i; 
  }

  delete [] scan; 
  return ret; 
} 

void Node::partition(Point *P, int *swap, int *bin, int *map, int threshold, 
                     double corner_x, double corner_y, double corner_z, 
                     double size) {
  int num_points = last_ - first_ + 1; 
  bool is_leaf = (num_points <= threshold ? true : false); 

  // Set complete_ LCO of the parent if the current node is not nullptr
  if (parent_) 
    hpx_call_when(complete_, parent_->complete(), hpx_lco_set_action, 
                  HPX_NULL, NULL, 0); 

  if (is_leaf) {
    // Set complete_ LCO of itself
    hpx_lco_and_set_num(complete_, 8, HPX_NULL); 
  } else {
    // Compute center of the node 
    double h = size / pow(2, idx_.level()); 
    double center_x = corner_x + (idx_.x() + 0.5) * h; 
    double center_y = corner_y + (idx_.y() + 0.5) * h; 
    double center_z = corner_z + (idx_.z() + 0.5) * h; 

    // Sort points into 8 octants of the node
    int stat[8] = {0}; 
    int offset[8] = {0}; 

    for (int i = first_; i <= last_; ++i) {
      Point *j = &P[map[i]]; 
      int temp = 4 * (j->z() > center_z) + 2 * (j->y() > center_y) + 
        (j->x() > center_x); 
      bin[i] = temp; 
      stat[temp]++;
    }

    offset[0] = first_; 
    for (int i = 1; i < 8; ++i)
      offset[i] = offset[i - 1] + stat[i - 1]; 

    // Rearrange map
    for (int i = first_; i <= last_; ++i) {
      int j = offset[bin[i]]++; 
      swap[j] = map[i]; 
    }

    for (int i = first_; i <= last_; ++i) 
      map[i] = swap[i]; 

    // Create child nodes
    for (int i = 0; i < 8; ++i) {
      if (stat[i]) {
        Node *child = new Node{idx_.child(i), offset[i] - stat[i], 
                               offset[i] - 1, this}; 
        child_[i] = child; 
        
        hpx_call(HPX_HERE, partition_node_action, HPX_NULL, 
                 &child, &P, &swap, &bin, &map); 
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
