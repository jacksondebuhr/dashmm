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


#if 0

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


#endif


/////////////////////////////////////////////////////////////////////
// Moving to a Node class
/////////////////////////////////////////////////////////////////////


hpx_action_t Node::partition_ = HPX_ACTION_NULL;


Node::~Node() {
  if (left) {
    delete left;
  }
  if (right) {
    delete right;
  }
}


void Node::partition(Particle *parts, int n_parts, int thresh,
                     hpx_addr_t sync) {
  hpx_addr_t done = sync;
  if (sync == HPX_NULL) {
    done = hpx_lco_future_new(0);
    assert(done != HPX_NULL);
  }

  Node *thisptr = this;
  hpx_call(HPX_HERE, partition_, done, &thisptr, &parts, &n_parts, &thresh);
}



/////////////////////////////////////////////////////////////////////
// Node actions
/////////////////////////////////////////////////////////////////////


void Node::register_actions() {
  HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                      Node::partition_,
                      Node::partition_handler,
                      HPX_POINTER, HPX_POINTER, HPX_INT, HPX_INT);
}


int Node::partition_handler(Node *node, Particle *parts, int n_parts,
                            int thresh) {
  if (n_parts <= thresh) {
    node->parts = parts;
    node->count = n_parts;
    return HPX_SUCCESS;
  } else {
    // find split point for the particles
    Particle *P_begin{parts};
    Particle *P_end{&P_begin[n_parts]};

    double x_split = 0.5 * (node->low + node->high);
    auto x_comp = [x_split](Particle &a) {
      return a.pos < x_split;
    };
    Particle *splitter = std::partition(P_begin, P_end, x_comp);
    int n_parts_left = splitter - P_begin;
    int n_parts_right = P_end - splitter;

    // Completion detection
    hpx_addr_t part_done = hpx_lco_and_new(2);
    assert(part_done != HPX_NULL);

    // call partition on the children, using that reduction as the result
    if (n_parts_left) {
      node->left = new Node(node->low, x_split);
      node->left->partition(P_begin, n_parts_left, thresh, part_done);
    } else {
      hpx_lco_and_set(part_done, HPX_NULL);
    }

    if (n_parts_right) {
      node->right = new Node(x_split, node->high);
      node->right->partition(splitter, n_parts_right, thresh, part_done);
    } else {
      hpx_lco_and_set(part_done, HPX_NULL);
    }

    // This passes completion information up the chain to the top, as well
    // as deleting the synchronization objects.
    return hpx_call_when_cc(part_done, part_done, hpx_lco_delete_action,
                            nullptr, 0);
  }
}

