#include <algorithm>

#include "tree.h"
#include "particle.h"


struct Node_data {
  hpx_addr_t left_;
  hpx_addr_t right_;
  double low_;
  double high_;
  moment_t moment_;
  hpx_addr_t parts_;
  int first_;
  int last_;
};


/////////////////////////////////////////////////////////////////////
// Private "Members"
/////////////////////////////////////////////////////////////////////


bool Node_is_leaf(Node_data *local) {
  return (local->left_ == HPX_NULL && local->right_ == HPX_NULL);
}


double Node_compute_direct(Node_data *local, double pos, int part) {
  Particle *P{nullptr};
  assert(hpx_gas_try_pin(local->parts_, (void **)&P));
  
  double retval{0.0};
  for (int i = local->first_; i != local->last_; ++i) {
    if (i != part) {
      double dist{fabs(pos - P[i].x())};
      retval -= P[i].m() / (dist + 1.0e-10);
    }
  }
  
  hpx_gas_unpin(local->parts_);
  return retval;
}


double Node_compute_approx(Node_data *local, double pos) {
  double dist{fabs(pos - local->moment_.xcom_)};
  double retval{-local->moment_.mtot_ / dist};
  retval -= local->moment_.Q00_ / (2.0 * dist * dist * dist);
  return retval;
}



/////////////////////////////////////////////////////////////////////
// Actions etc...
/////////////////////////////////////////////////////////////////////


int Node_construction(Node_data *local, double low, double high) {
  local->low_ = low;
  local->high_ = high;
  local->left_ = HPX_NULL;
  local->right_ = HPX_NULL;
  local->parts_ = HPX_NULL;
  local->first_ = 0;
  local->last_ = 0;
  local->moment_ = moment_t{0.0, 0.0, 0.0};
  
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_PINNED, 
           Node_construction_action, Node_construction,
           HPX_POINTER, HPX_DOUBLE, HPX_DOUBLE);


//Gee, isn't this stupid. I need to call the action from the action so I 
// have to do this nonsense. Even if I did the registration myself, and I
// declared the action, I would still have to have the prototype ready, and
// so that would not actually save me this repetition. Awesome!
int Node_destruction(Node_data *local, void *UNUSED, size_t UNWANTED);
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED | HPX_PINNED,
           Node_destruction_action, Node_destruction,
           HPX_POINTER, HPX_POINTER, HPX_SIZE_T);

int Node_destruction(Node_data *local, void *UNUSED, size_t UNWANTED) {
  hpx_addr_t done{hpx_lco_and_new(2)};
  assert(done != HPX_NULL);
  
  if (local->left_ != HPX_NULL) {
    hpx_call(local->left_, Node_destruction_action, done, nullptr, 0);
  } else {
    hpx_lco_set(done, 0, nullptr, HPX_NULL, HPX_NULL);
  }
  if (local->right_ != HPX_NULL) {
    hpx_call(local->right_, Node_destruction_action, done, nullptr, 0);
  } else {
    hpx_lco_set(done, 0, nullptr, HPX_NULL, HPX_NULL);
  }
  
  hpx_lco_wait(done);
  hpx_lco_delete(done, HPX_NULL);
  
  return HPX_SUCCESS;
}


//This seems to be needed because the variadic calls and things somehow
// are not able to take no arguments, for things like hpx_lco_delete_action.
int thatsrightthisisstupid(void *DUMB, size_t DUMBER) {
  hpx_addr_t target{hpx_thread_current_target()};
  hpx_lco_delete(target, HPX_NULL);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
           thatsrightthisisstupidaction, thatsrightthisisstupid,
           HPX_POINTER, HPX_SIZE_T);


int Node_partition(Node_data *local, hpx_addr_t parts, 
                            int first, int last, int limit) {
  local->first_ = first;
  local->last_ = last;
  local->parts_ = parts;
  
  if ((last - first) <= limit) {
    return HPX_SUCCESS;
  }
  
  double span = 0.5 * (local->high_ - local->low_);
  double divide = local->low_ + span;
  
  Node left(local->low_, divide);
  Node right(divide, local->high_);
  left.drop_ownership();
  right.drop_ownership();
  local->left_ = left.data();
  local->right_ = right.data();
  
  //find separator
  Particle *P{nullptr};
  assert(hpx_gas_try_pin(parts, (void **)&P));
  Particle comparator(local->low_ + span, 0.0);
  Particle *separator = std::lower_bound(&P[first], &P[last], comparator);
  int split = separator - &P[first] + first;
  hpx_gas_unpin(parts);
  
  //recurse
  hpx_addr_t part_done = hpx_lco_and_new(2);
  assert(part_done != HPX_NULL);
  left.partition(parts, first, split, limit, part_done);
  right.partition(parts, split, last, limit, part_done);
  
  //then do a call when cc on the part_done LCO to delete it
  hpx_call_when_cc(part_done, part_done, thatsrightthisisstupidaction, 
                   nullptr, nullptr, nullptr, 0);
}
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           Node_partition_action, Node_partition,
           HPX_POINTER, HPX_ADDR, HPX_INT, HPX_INT, HPX_INT);


void moment_reduction_ident_handler(moment_t *ident, size_t bytes) {
  ident->mtot_ = 0.0;
  ident->xcom_ = 0.0;
  ident->Q00_  = 0.0;
}
HPX_ACTION(HPX_FUNCTION, 0, 
           moment_reduction_ident,
           moment_reduction_ident_handler, HPX_POINTER, HPX_SIZE_T);


void moment_reduction_op_handler(moment_t *lhs, 
                                 const moment_t *rhs, size_t bytes) {
  double newtot = lhs->mtot_ + rhs->mtot_;
  double newcom{0.0};
  double newQ00{0.0};
  if (lhs->mtot_ > 0.0) {
    newcom = lhs->mtot_ * lhs->xcom_ + rhs->mtot_ * rhs->xcom_;
    newcom /= newtot;
    
    double dxl{lhs->xcom_ - newcom};
    newQ00 = 2.0 * lhs->mtot_ * dxl * dxl + lhs->Q00_;
    
    double dxr{rhs->xcom_ - newcom};
    newQ00 += 2.0 * rhs->mtot_ * dxr * dxr + rhs->Q00_;
  }
  
  lhs->mtot_ = newtot;
  lhs->xcom_ = newcom;
  lhs->Q00_  = newQ00;
}
HPX_ACTION(HPX_FUNCTION, 0, 
           moment_reduction_op,
           moment_reduction_op_handler, HPX_POINTER, HPX_POINTER, HPX_SIZE_T);


int Node_moment_reduction_finished(void *UNUSED, size_t UNWANTED) {
  hpx_addr_t target = hpx_thread_current_target();
  moment_t contval{};
  hpx_lco_get(target, sizeof(contval), &contval);
  hpx_lco_delete(target, HPX_NULL);
  HPX_THREAD_CONTINUE(contval);
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED, 
           Node_moment_reduction_finished_action,
           Node_moment_reduction_finished,
           HPX_POINTER, HPX_SIZE_T);


int Node_compute_moment(Node_data *local, void *UNUSED, size_t UNWANTED) {
  if (Node_is_leaf(local)) {
    //This is S->M
    Particle *P{nullptr};
    assert(hpx_gas_try_pin(local->parts_, (void **)&P));
    
    for (int i = local->first_; i != local->last_; ++i) {
      local->moment_.mtot_ += P[i].m();
      local->moment_.xcom_ += P[i].m() * P[i].x();
    }
    if (local->moment_.mtot_ > 0.0) {
      local->moment_.xcom_ /= local->moment_.mtot_;
      for (int i = local->first_; i != local->last_; ++i) {
        double dx = P[i].x() - local->moment_.xcom_;
        local->moment_.Q00_ += 2.0 * P[i].m() * dx * dx;
      }
    }
    
    hpx_gas_unpin(local->parts_);
    
    //done, continue the moment value
    HPX_THREAD_CONTINUE(local->moment_);
  } else {
    Node left(local->left_);
    Node right(local->right_);
    
    hpx_addr_t reduce_moment = hpx_lco_reduce_new(2, sizeof(moment_t), 
                                    moment_reduction_ident,
                                    moment_reduction_op);
    assert(reduce_moment != HPX_NULL);
    
    left.compute_moment(reduce_moment);
    right.compute_moment(reduce_moment);
    
    //set up, so do a call when cc
    hpx_call_when_cc(reduce_moment, reduce_moment, 
                     Node_moment_reduction_finished_action, 
                     nullptr, nullptr, nullptr, 0);
  }
}
HPX_ACTION(HPX_DEFAULT, HPX_PINNED | HPX_MARSHALLED,
           Node_compute_moment_action, Node_compute_moment,
           HPX_POINTER, HPX_POINTER, HPX_SIZE_T);


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


int Node_potential_finished(hpx_addr_t whendone) {
  hpx_addr_t target = hpx_thread_current_target();
  double contval{0.0};
  hpx_lco_get(target, sizeof(contval), &contval);
  hpx_lco_delete(target, HPX_NULL);
  if (whendone == HPX_NULL) {
    HPX_THREAD_CONTINUE(contval);
  } else {
    //Only the root will (possibly) have a nonzero whendone
    hpx_thread_continue(&contval, &whendone);
  }
}
HPX_ACTION(HPX_DEFAULT, 0,
           Node_potential_finished_action, Node_potential_finished,
           HPX_ADDR);


int Node_compute_potential(Node_data *local, double pos, 
                           int part, double theta_c, hpx_addr_t whendone) {
  double dist{fabs(pos - local->moment_.xcom_)};
  double size{local->high_ - local->low_};
  double theta{size / dist};
  if (theta < theta_c) {
    double retval{Node_compute_approx(local, pos)};
    HPX_THREAD_CONTINUE(retval);
  } else {
    if (Node_is_leaf(local)) {
      double retval{Node_compute_direct(local, pos, part)};
      HPX_THREAD_CONTINUE(retval);
    } else {
      //spawn work at children, and set up call whens to make the results
      // get continued
      Node left(local->left_);
      Node right(local->right_);
      
      hpx_addr_t potential = hpx_lco_reduce_new(2, sizeof(double), 
                                double_sum_ident, double_sum_op);
      assert(potential != HPX_NULL);
      
      left.compute_potential(pos, part, theta_c, potential);
      right.compute_potential(pos, part, theta_c, potential);
      
      hpx_call_when_cc(potential, potential, Node_potential_finished_action,
                       nullptr, nullptr, &whendone);
    }
  }
}
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           Node_compute_potential_action, Node_compute_potential,
           HPX_POINTER, HPX_DOUBLE, HPX_INT, HPX_DOUBLE, HPX_ADDR);



/////////////////////////////////////////////////////////////////////
// Interface to Node
/////////////////////////////////////////////////////////////////////

Node::Node(hpx_addr_t data) {
  data_ = data;
  owner_ = false;
}


Node::Node(double lowbound, double highbound, hpx_addr_t colocate) {
  hpx_addr_t where {colocate};
  if (where == HPX_NULL) {
    where = HPX_HERE;
  }
  data_ = hpx_gas_alloc_local_at_sync(sizeof(Node_data), 0, where);
  
  if (data_ != HPX_NULL) {
    hpx_call_sync(data_, Node_construction_action, nullptr, 0,
                    &lowbound, &highbound);
  }
}


Node::Node(const Node &other) {
  data_ = other.data_;
  owner_ = false;
}


Node::Node(Node &&other) {
  data_ = other.data_;
  owner_ = other.owner_;
  other.owner_ = false;
}


Node::~Node() {
  if (owner_) {
    hpx_call_sync(data_, Node_destruction_action, nullptr, 0, nullptr, 0);
  
    hpx_gas_free_sync(data_);
  }
}


Node &Node::operator=(Node &&other) {
  if (owner_) {
    hpx_gas_free(data_, HPX_NULL);
  }
  data_ = other.data_;
  owner_ = other.owner_;
  other.owner_ = false;
  
  return *this;
}


Node &Node::operator=(const Node &other) {
  data_ = other.data_;
  owner_ = false;
  return *this;
}

  
void Node::partition_sync(hpx_addr_t parts, int first, int last, int limit) {
  hpx_addr_t done = hpx_lco_future_new(0);
  assert(done != HPX_NULL);
  partition(parts, first, last, limit, done);
  hpx_lco_wait(done);
  hpx_lco_delete(done, HPX_NULL);
}


hpx_addr_t Node::partition(hpx_addr_t parts, int first, int last, int limit, 
                       hpx_addr_t sync) {
  hpx_call(data_, Node_partition_action, sync, &parts, &first, &last, &limit);
  return sync;
}

  
moment_t Node::compute_moment_sync() {
  hpx_addr_t done = hpx_lco_future_new(sizeof(moment_t));
  assert(done != HPX_NULL);
  compute_moment(done);
  moment_t retval{};
  hpx_lco_get(done, sizeof(retval), &retval);
  hpx_lco_delete(done, HPX_NULL);
  return retval;
}


hpx_addr_t Node::compute_moment(hpx_addr_t sync) {
  hpx_call(data_, Node_compute_moment_action, sync, nullptr, 0);
  return sync;
}
 
  
double Node::compute_potential_sync(double pos, int part, double theta_c) {
  hpx_addr_t done = hpx_lco_future_new(sizeof(double));
  assert(done != HPX_NULL);
  compute_potential(pos, part, theta_c, done);
  double retval{0.0};
  hpx_lco_get(done, sizeof(retval), &retval);
  hpx_lco_delete(done, HPX_NULL);
  return retval;
}


hpx_addr_t Node::compute_potential(double pos, int part, double theta_c, 
                               hpx_addr_t sync) {
  hpx_addr_t empty{HPX_NULL};
  hpx_call(data_, Node_compute_potential_action, sync, 
                                &pos, &part, &theta_c, &empty);
  return sync;
}


void Node::compute_potential_and_continue(double pos, int part, double theta_c,
                                      SetApproxContinuation cont,
                                      hpx_addr_t sync) {
  assert(sync != HPX_NULL);
  hpx_call_with_continuation(data_, Node_compute_potential_action, 
                         cont.target(), cont.action(),
                         &pos, &part, &theta_c, &sync);
}


