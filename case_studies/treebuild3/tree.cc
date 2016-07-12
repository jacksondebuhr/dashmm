#include <cstdlib>
#include <cmath>
#include <cstdint> 
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
  hpx_bcast_rsync(set_points_action, &src, &tar, &nsrc, &ntar, &datatype); \
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

  // Choose a coarser level L to perform a uniform partition, where 
  // 8^L >= num_ranks
  int L = ceil(log(num_ranks) / log(8)); 
  int num_grids = pow(8, L); 

  hpx_addr_t coarse_count = hpx_gas_calloc_cyclic(num_ranks, 
                                                  sizeof(hpx_addr_t), 0); 
  hpx_addr_t unif_grid = hpx_gas_calloc_cyclic(num_ranks, 
                                               sizeof(Node) * num_grids * 2, 
                                               0); 
  hpx_addr_t sorted_src = hpx_gas_calloc_cyclic(num_ranks, sizeof(Point *), 0);
  hpx_addr_t sorted_tar = hpx_gas_calloc_cyclic(num_ranks, sizeof(Point *), 0); 

  // Broadcast addresses to perform coarse level partition 
  hpx_bcast_rsync(init_coarse_partition_action, &coarse_count, &unif_grid, &L);


  // Sort the source and target points on each rank into coarse grids
  hpx_bcast_rsync(coarse_partition_action, &src, &nsrc, &tar, &ntar, 
                  &coarse_count, &unif_grid, &sorted_src, &sorted_tar, 
                  &L, &corner_x, &corner_y, &corner_z, &size); 







  /*



  // Delete trees 
  hpx_bcast_rsync(delete_tree_action, &src_tree_root, &tar_tree_root); 

  hpx_lco_delete_sync(domain_geometry); 
  hpx_gas_free_sync(src_tree_root); 
  hpx_gas_free_sync(tar_tree_root);
  */
  hpx_exit(HPX_SUCCESS); 
}
HPX_ACTION(HPX_DEFAULT, 0, partition_points_action, partition_points_handler, 
           HPX_INT, HPX_ADDR, HPX_INT, HPX_ADDR, HPX_INT); 


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

int init_coarse_partition_handler(hpx_addr_t coarse_count, 
                                  hpx_addr_t unif_grid, int level) {
  int rank = hpx_get_my_rank(); 
  int num_ranks = hpx_get_num_ranks(); 
  int num_per_dim = pow(2, level); 
  int num_grids = pow(8, level); 

  // Setup user LCO count 
  hpx_addr_t curr_coarse_count 
    = hpx_addr_add(count, sizeof(hpx_addr_t) * rank, sizeof(hpx_addr_t)); 
  hpx_addr_t *user_lco{nullptr}; 
  if (!hpx_gas_try_pin(curr_coarse_count, (void **)&user_lco))
    return HPX_ERROR; 
  
  int init[2] = {num_ranks, num_grids * 2}; 
  *user_lco = hpx_lco_user_new(sizeof(int) * (2 + num_grids * 2), 
                               coarse_grid_count_init_action, 
                               coarse_grid_count_op_action, 
                               coarse_grid_count_predicate_action,
                               &init, sizeof(init)); 
  hpx_gas_unpin(curr_coarse_count); 
 
  // Setup the uniform grid 
  hpx_addr_t curr_unif_grid 
    = hpx_addr_add(unif_grid, sizeof(Node) * num_grids * 2 * rank, 
                   sizeof(Node) * num_grids * 2); 
  Node *n{nullptr}, *s{nullptr}, *t{nullptr}; 

  if (!hpx_gas_try_pin(curr_unif_grid, (void **)&n))
    return HPX_ERROR; 

  s = &n[0]; 
  t = &n[num_grids]; 

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
  hpx_gas_unpin(curr_unif_grid); 

  return HPX_SUCCESS; 
}
HPX_ACTION(DEFAULT, 0, init_partition_action, init_partition_handler, 
           HPX_ADDR, HPX_ADDR, HPX_INT); 

int coarse_partition_handler(hpx_addr_t src, int nsrc, hpx_addr_t tar, int ntar, 
                             hpx_addr_t coarse_count, hpx_addr_t unif_grid, 
                             hpx_addr_t sorted_src, hpx_addr_t sorted_tar, 
                             int L, double corner_x, double corner_y, 
                             double corner_z, double size) {
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
  
  int dim = pow(2, L), dim2 = dim * dim, dim3 = dim2 * dim; 
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

  hpx_gas_unpin(curr_src); 
  hpx_gas_unpin(curr_tar); 

  // Exchange counts to determine a distribution of the coarse level 
  for (int i = 0; i < num_ranks; ++i) {
    hpx_parcel_t *p = hpx_parcel_acquire(local_cnt, sizeof(int) * dim3 * 2);
    hpx_addr_t count_i = hpx_addr_add(coarse_count, sizeof(hpx_addr_t) * i, 
                                      sizeof(hpx_addr_t)); 
    hpx_parcel_set_target(p, count_i); 
    hpx_parcel_set_action(p, exchange_coarse_grid_count_action); 
    hpx_parcel_send(p, HPX_NULL); 
  }

  hpx_addr_t curr_count = hpx_addr_add(coarse_count, sizeof(hpx_addr_t) * rank, 
                                       sizeof(hpx_addr_t)); 
  void *user_lco_buffer{nullptr}; 
  int *global_cnt{nullptr}; 
  hpx_lco_getref(curr_count, 1, &user_lco_buffer); 
  global_cnt = reinterpret_cast<int *>(static_cast<char *>(user_lco_buffer) 
                                       + sizeof(int) * 2); 

  // Compute point distribution 
  int *distribution = distribute_points(num_ranks, global_cnt, dim3); 
 
  // Prepare to exchange points 
  hpx_addr_t curr_sorted_src = hpx_addr_add(sorted_src, sizeof(Point *) * rank, 
                                            sizeof(Point *)); 
  hpx_addr_t curr_sorted_tar = hpx_addr_add(sorted_tar, sizeof(Point *) * rank, 
                                            sizeof(Point *)); 
  Point **sorted_s{nullptr}, **sorted_t{nullptr}; 

  if (!hpx_gas_try_pin(curr_sorted_src, (void **)&sorted_s) || 
      !hpx_gas_try_pin(curr_sorted_tar, (void **)&sorted_t)) 
    return HPX_ERROR; 

  





hpx_addr_t unif_grid, 
                             hpx_addr_t sorted_src, hpx_addr_t sorted_tar, 



  // Compute number of points owned by the current rank
  int start = (rank == 0 ? 0 : distribution[rank - 1] + 1); 
  int end = distribution[rank]; 
  int sorted_nsrc{0}, sorted_ntar{0}; 

  for (int i = start; i <= end; ++i) {
    sorted_nsrc += global_cnt[i]; 
    sorted_ntar += global_cnt[dim3 + i]; 
  }



  
  
  
  
  return HPX_SUCCESS; 
}
HPX_ACTION(HPX_DEFAULT, 0, coarser_partition_action, coarser_partition_handler, 
           HPX_ADDR, HPX_INT, HPX_ADDR, HPX_INT, HPX_ADDR, HPX_ADDR, 
           HPX_ADDR, HPX_ADDR, HPX_INT, HPX_DOUBLE, HPX_DOUBLE, 
           HPX_DOUBLE, HPX_DOUBLE); 

int exchange_coarse_grid_count_handler(void *local_target, 
                                       void *args, size_t size) {
  // local_target is (hpx_addr_t *) type pointing to user-defined LCO 
  hpx_addr_t user_lco = 
    *(reinterpret_cast<hpx_addr_t *>(static_cast<char *>(local_target))); 
  hpx_lco_set(user_lco, size, args, HPX_NULL, HPX_NULL); 
  return HPX_SUCCESS; 
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED, exchange_coarse_grid_count_action, 
           exchange_coarse_grid_count_handler, HPX_POINTER, HPX_SIZE_T); 


void coarse_grid_count_init_handler(void *data, const size_t size, 
                                    int *init, size_t init_size) {

  int *input = reinterpret_cast<int *>(static_cast<char *>(data)); 
  
  input[0] = init[0]; // # of inputs; 
  input[1] = init[1]; // # of terms 

  for (int i = 0; i < init[1]; ++i) 
    input[i + 2] = 0; 
} 
HPX_ACTION(HPX_FUNCTION, 0, coarse_grid_count_init_action, 
           coarse_grid_count_init_handler, HPX_POINTER, HPX_SIZE_T, 
           HPX_POINTER, HPX_SIZE_T); 

void coarse_grid_count_op_handler(void *data, int *rhs, size_t size) {
  int *input = reinterpret_cast<int *>(static_cast<char *>(data)); 
  int *lhs = &input[2]; 

  input[0]--; 
  int nterms = input[1]; 
  for (int i = 0; i < nterms; ++i) 
    lhs[i] += rhs[i]; 
}
HPX_ACTION(HPX_FUNCTION, 0, coarse_grid_count_op_action, 
           coarse_grid_count_op_handler, HPX_POINTER, HPX_POINTER, HPX_SIZE_T); 

bool coarse_grid_count_predicate_handler(void *data, size_t size) {
  int *input = reinterpret_cast<int *>(static_cast<char *>(data)); 
  return (input[0] == 0); 
}
HPX_ACTION(HPX_FUNCTION, 0, coarse_grid_count_predicate_action, 
           coarse_grid_count_predicate_handler, HPX_POINTER, HPX_SIZE_T); 

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

