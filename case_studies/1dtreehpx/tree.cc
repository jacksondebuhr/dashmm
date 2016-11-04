// ============================================================================
//  High Performance ParalleX Library (hpx-apps)
//
//  Copyright (c) 2013-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license.  See the COPYING file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// ============================================================================

#include "tree.h"

#include <algorithm>
#include <cmath>

#include "hpx/hpx.h"

#include "particle.h"


static double domain_size_;
static int domain_count_;


struct Moment {
  double mtot;
  double xcom;
  double Q00;
};


struct Node {
  hpx_addr_t left;
  hpx_addr_t right;
  double low;
  double high;
  Moment moments;
  hpx_addr_t parts;
  int count;
};


/////////////////////////////////////////////////////////////////////
// Internal Utilities
/////////////////////////////////////////////////////////////////////


double domain_low_bound(int which) {
  return (domain_size_ / domain_count_) * which;
}


double domain_high_bound(int which) {
  return (domain_size_ / domain_count_) * (which + 1);
}


int map_bounds_to_locality(double low, double high) {
  double domain_span = domain_size_ / domain_count_;
  int low_idx = (int)floor(low / domain_span);
  int high_idx = (int)floor(high / domain_span);
  assert(high_idx >= low_idx);

  //Easiest case first - only one choice really
  if (low_idx == high_idx) {
    return low_idx;
  }

  //harder case - pick the one with the most overlap
  if (high_idx - low_idx == 1) {
    double delta_low = domain_low_bound(high_idx) - low;
    double delta_high = high - domain_low_bound(high_idx);
    return (delta_high > delta_low ? high_idx : low_idx);
  }

  //hardest case - that means the most possible answers. We just pick the
  // first. There is a potential for a better choice here.
  return low_idx;
}


void compute_moments(Node *node, Particle *P, int n_parts) {
  node->moments.mtot = 0.0;
  node->moments.xcom = 0.0;
  node->moments.Q00 = 0.0;
  for (int i = 0; i < n_parts; ++i) {
    node->moments.mtot += P[i].mass;
    node->moments.xcom += P[i].pos * P[i].mass;
  }
  if (node->moments.mtot > 0.0) {
    node->moments.xcom /= node->moments.mtot;
    for (int i = 0; i < n_parts; ++i) {
      double dx = P[i].pos - node->moments.xcom;
      node->moments.Q00 += 2.0 * P[i].mass * dx * dx;
    }
  }
}


bool is_leaf(Node *node) {
  return node->left == HPX_NULL && node->right == HPX_NULL;
}


double node_compute_approx(Node *node, double pos) {
  double dist = fabs(pos - node->moments.xcom);
  double retval = -node->moments.mtot / dist;
  retval -= node->moments.Q00 / (2.0 * dist * dist * dist);
  return retval;
}


double node_compute_direct(Particle *P, int count, double pos) {
  double retval{0.0};
  for (int i = 0; i < count; ++i) {
    double dist = fabs(pos - P[i].pos);
    if (dist > 0) {
      retval -= P[i].mass / dist;
    }
  }
  return retval;
}


/////////////////////////////////////////////////////////////////////
// Actions
/////////////////////////////////////////////////////////////////////


int broadcast_domain_handler(double size) {
  domain_size_ = size;
  domain_count_ = hpx_get_num_ranks();
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           broadcast_domain_action, broadcast_domain_handler,
           HPX_DOUBLE);


void moment_reduction_ident_handler(Moment *ident, size_t bytes) {
  ident->mtot = 0.0;
  ident->xcom = 0.0;
  ident->Q00  = 0.0;
}
HPX_ACTION(HPX_FUNCTION, 0,
           moment_reduction_ident, moment_reduction_ident_handler,
           HPX_POINTER, HPX_SIZE_T);


void moment_reduction_op_handler(Moment *lhs,
                                 const Moment *rhs, size_t bytes) {
  double newtot = lhs->mtot + rhs->mtot;
  double newcom{0.0};
  double newQ00{0.0};
  if (newtot > 0.0) {
    newcom = lhs->mtot * lhs->xcom + rhs->mtot * rhs->xcom;
    newcom /= newtot;

    double dxl{lhs->xcom - newcom};
    newQ00 = 2.0 * lhs->mtot * dxl * dxl + lhs->Q00;

    double dxr{rhs->xcom - newcom};
    newQ00 += 2.0 * rhs->mtot * dxr * dxr + rhs->Q00;
  }

  lhs->mtot = newtot;
  lhs->xcom = newcom;
  lhs->Q00  = newQ00;
}
HPX_ACTION(HPX_FUNCTION, 0,
           moment_reduction_op, moment_reduction_op_handler,
           HPX_POINTER, HPX_POINTER, HPX_SIZE_T);


int save_and_continue_moments_handler(Node *node, hpx_addr_t moments) {
  hpx_lco_get(moments, sizeof(Moment), &node->moments);
  hpx_lco_delete_sync(moments);
  return HPX_THREAD_CONTINUE(node->moments);
}
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           save_and_continue_moments_action, save_and_continue_moments_handler,
           HPX_POINTER, HPX_ADDR);


int partition_node_handler(Node *node, hpx_addr_t parts, int n_parts,
                           int n_partition);
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           partition_node_action, partition_node_handler,
           HPX_POINTER, HPX_ADDR, HPX_INT, HPX_INT);

int partition_node_handler(Node *node, hpx_addr_t parts, int n_parts,
                           int n_partition) {
  if (n_parts <= n_partition) {
    node->parts = parts;
    node->count = n_parts;

    Particle *P;
    assert(hpx_gas_try_pin(parts, (void **)&P));
    compute_moments(node, P, n_parts);
    hpx_gas_unpin(parts);

    return HPX_THREAD_CONTINUE(node->moments);
  } else {
    // find split point for the particles particles
    Particle *P_begin;
    assert(hpx_gas_try_pin(parts, (void **)&P_begin));
    Particle *P_end = &P_begin[n_parts];

    double x_split = 0.5 * (node->low + node->high);
    auto x_comp = [x_split](Particle &a) {
      return a.pos < x_split;
    };
    Particle *splitter = std::partition(P_begin, P_end, x_comp);
    int n_parts_left = splitter - P_begin;
    int n_parts_right = P_end - splitter;
    hpx_gas_unpin(parts);

    //get a reduction to combine moments once partitioning is finished
    hpx_addr_t reduce = hpx_lco_reduce_new(2, sizeof(Moment),
                          moment_reduction_ident, moment_reduction_op);
    assert(reduce != HPX_NULL);

    // call partition on the children, using that reduction as the result
    if (n_parts_left) {
      node->left = create_node(node->low, x_split);
      hpx_addr_t parts_left = hpx_gas_alloc_local_at_sync(1,
                    sizeof(Particle) * n_parts_left, 0, node->left);
      hpx_addr_t cpy_done = hpx_lco_future_new(0);
      assert(cpy_done != HPX_NULL);
      hpx_gas_memcpy(parts_left, parts, sizeof(Particle) * n_parts_left,
                     cpy_done);
      hpx_lco_wait(cpy_done);
      hpx_lco_delete_sync(cpy_done);
      hpx_call(node->left, partition_node_action, reduce,
               &parts_left, &n_parts_left, &n_partition);
    } else {
      Moment empty{0.0, 0.0, 0.0};
      hpx_lco_set_rsync(reduce, sizeof(Moment), &empty);
    }

    if (n_parts_right) {
      node->right = create_node(x_split, node->high);
      hpx_addr_t parts_right = hpx_gas_alloc_local_at_sync(1,
                    sizeof(Particle) * n_parts_right, 0, node->right);
      hpx_addr_t from_addx = hpx_addr_add(parts,
                                          sizeof(Particle) * n_parts_left,
                                          sizeof(Particle) * n_parts);
      hpx_addr_t cpy_done = hpx_lco_future_new(0);
      assert(cpy_done != HPX_NULL);
      hpx_gas_memcpy(parts_right, from_addx,
                     sizeof(Particle) * n_parts_right, cpy_done);
      hpx_lco_wait(cpy_done);
      hpx_lco_delete_sync(cpy_done);
      hpx_call(node->right, partition_node_action, reduce,
               &parts_right, &n_parts_right, &n_partition);
    } else {
      Moment empty{0.0, 0.0, 0.0};
      hpx_lco_set_rsync(reduce, sizeof(Moment), &empty);
    }

    return hpx_call_when_cc(reduce, hpx_thread_current_target(),
                            save_and_continue_moments_action,
                            &reduce);
  }
}


/////////////////////////////////////////////////////////////////////////////


void double_sum_ident_handler(double *id, size_t UNUSED) {
  *id = 0.0;
}
HPX_ACTION(HPX_FUNCTION, 0,
           double_sum_ident, double_sum_ident_handler, HPX_POINTER, HPX_SIZE_T);


void double_sum_op_handler(double *lhs, const double *rhs, size_t UNUSED) {
  *lhs += *rhs;
}
HPX_ACTION(HPX_FUNCTION, 0,
           double_sum_op, double_sum_op_handler,
           HPX_POINTER, HPX_POINTER, HPX_SIZE_T);


int particle_set_approx_handler(Particle *part, hpx_addr_t phi) {
  hpx_lco_get(phi, sizeof(double), &part->phi);
  hpx_lco_delete_sync(phi);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           particle_set_approx_action, particle_set_approx_handler,
           HPX_POINTER, HPX_ADDR);


int potential_reduction_complete_handler(void *UNUSED, size_t bytes) {
  hpx_addr_t target = hpx_thread_current_target();
  double contval{0.0};
  hpx_lco_get(target, sizeof(contval), &contval);
  hpx_lco_delete_sync(target);
  return HPX_THREAD_CONTINUE(contval);
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
           potential_reduction_complete_action,
           potential_reduction_complete_handler,
           HPX_POINTER, HPX_SIZE_T);


int node_compute_potential_handler(Node *node, double pos, double theta);
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           node_compute_potential_action, node_compute_potential_handler,
           HPX_POINTER, HPX_DOUBLE, HPX_DOUBLE);

int node_compute_potential_handler(Node *node, double pos, double theta) {
  double dist = fabs(pos - node->moments.xcom);
  double size = node->high - node->low;
  double test = size / dist;
  if (test < theta) {
    double retval = node_compute_approx(node, pos);
    return HPX_THREAD_CONTINUE(retval);
  } else if (is_leaf(node)) {
    Particle *P{nullptr};
    assert(hpx_gas_try_pin(node->parts, (void **)&P));
    double retval = node_compute_direct(P, node->count, pos);
    hpx_gas_unpin(node->parts);
    return HPX_THREAD_CONTINUE(retval);
  } else {
    //A reduction sums the results from the children
    hpx_addr_t potential = hpx_lco_reduce_new(2, sizeof(double),
                                  double_sum_ident, double_sum_op);
    assert(potential != HPX_NULL);

    //call the work on the children
    if (node->left != HPX_NULL) {
      hpx_call(node->left, node_compute_potential_action, potential,
                                                &pos, &theta);
    } else {
      double ident{0.0};
      hpx_lco_set(potential, sizeof(double), &ident, HPX_NULL, HPX_NULL);
    }
    if (node->right != HPX_NULL) {
      hpx_call(node->right, node_compute_potential_action, potential,
                                                &pos, &theta);
    } else {
      double ident{0.0};
      hpx_lco_set(potential, sizeof(double), &ident, HPX_NULL, HPX_NULL);
    }

    //set up work to do when this is done
    return hpx_call_when_cc(potential, potential,
                            potential_reduction_complete_action, nullptr, 0);
  }
}



int compute_and_save_handler(hpx_addr_t root, hpx_addr_t sync,
                             hpx_addr_t partaddx, double pos, double theta) {
  hpx_addr_t compdone = hpx_lco_future_new(sizeof(double));
  assert(compdone != HPX_NULL);
  hpx_call(root, node_compute_potential_action, compdone, &pos, &theta);
  hpx_call_when(compdone, partaddx, particle_set_approx_action, sync,
                                    &compdone);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           compute_and_save_action, compute_and_save_handler,
           HPX_ADDR, HPX_ADDR, HPX_ADDR, HPX_DOUBLE, HPX_DOUBLE);


int spawn_computation_handler(Node *node, hpx_addr_t root, hpx_addr_t sync,
                              double theta);
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           spawn_computation_action, spawn_computation_handler,
           HPX_POINTER, HPX_ADDR, HPX_ADDR, HPX_DOUBLE);

int spawn_computation_handler(Node *node, hpx_addr_t root, hpx_addr_t sync,
                              double theta) {
  if (is_leaf(node)) {
    Particle *P{nullptr};
    assert(hpx_gas_try_pin(node->parts, (void **)&P));

    for (int i = 0; i < node->count; ++i) {
      hpx_addr_t paddx = hpx_addr_add(node->parts, sizeof(Particle) * i,
                                      sizeof(Particle) * node->count);
      hpx_call(HPX_HERE, compute_and_save_action, HPX_NULL, &root, &sync,
                          &paddx, &P[i].pos, &theta);
    }

    hpx_gas_unpin(node->parts);
  } else {
    //not a leaf, so we continue the tree spawn
    if (node->left != HPX_NULL) {
      hpx_call(node->left, spawn_computation_action, HPX_NULL,
                           &root, &sync, &theta);
    }
    if (node->right != HPX_NULL) {
      hpx_call(node->right, spawn_computation_action, HPX_NULL,
                            &root, &sync, &theta);
    }
  }
  return HPX_SUCCESS;
}


/////////////////////////////////////////////////////////////////////
// Interface
/////////////////////////////////////////////////////////////////////


void broadcast_domain_size(double domain_size) {
  hpx_bcast_rsync(broadcast_domain_action, &domain_size);
}


hpx_addr_t create_node(double low, double high) {
  int where = map_bounds_to_locality(low, high);
  hpx_addr_t retval = hpx_gas_alloc_local_at_sync(1,
                          sizeof(Node), 0, HPX_THERE(where));
  Node vals{HPX_NULL, HPX_NULL, low, high,
              Moment{0.0, 0.0, 0.0}, HPX_NULL, 0};
  hpx_gas_memput_rsync(retval, &vals, sizeof(Node));
  return retval;
}


void partition_node_sync(hpx_addr_t node, hpx_addr_t parts, int n_parts,
                         int n_partition) {
  hpx_addr_t done = hpx_lco_future_new(sizeof(Moment));
  assert(done != HPX_NULL);
  hpx_call(node, partition_node_action, done, &parts, &n_parts, &n_partition);
  hpx_lco_wait(done);
  hpx_lco_delete_sync(done);
}


void spawn_potential_computation(hpx_addr_t root, hpx_addr_t sync,
                                 double theta) {
  hpx_call(root, spawn_computation_action, HPX_NULL, &root, &sync, &theta);
}
