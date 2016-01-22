// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license. See the LICENSE file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================


/// \file src/expansionref.cc
/// \brief Implementation of Expansion reference object


#include "include/expansionref.h"

#include <cassert>
#include <cstring>

#include <hpx/hpx.h>


namespace dashmm {


/////////////////////////////////////////////////////////////////////
// Stuff for the user LCO
/////////////////////////////////////////////////////////////////////


/// Part of the internal representation of the Expansion LCO
///
/// Expansion are user-defined LCOs. The data they contain are this object
/// and the serialized expansion. This object gives the number of expected
/// inputs, the number that have actually occurred, and a flag to indicate
/// if all of the expected inputs have been scheduled.
struct ExpansionLCOHeader {
  int arrived;
  int scheduled;
  int finished;
  size_t payload_size;
  char payload[];
};

/// Behavior codes for the Expansion LCO
///
/// The set operation for the Expansion LCO takes three forms. Two are simple:
/// incrementing the number count of scheduled inputs, and finalizing the
/// inputs. The third is for actually making contributions to the Expansion
/// data.
enum ExpansionLCOSetCodes {
  kFinish = 1,
  kSchedule = 2,
  kContribute = 3
};


/// Initialization handler for Expansion LCOs
///
/// This initialized an Expansion LCO given an input serialized
/// expansion. Often, this will just be the default constructed expansion,
/// but might be otherwise in specific cases.
void expansion_lco_init_handler(ExpansionLCOHeader *head, size_t bytes,
                                void *init, size_t init_bytes) {
  assert(bytes == init_bytes + sizeof(ExpansionLCOHeader));
  head->arrived = 0;
  head->scheduled = 0;
  head->finished = 0;
  head->payload_size = init_bytes;
  memcpy(head->payload, init, init_bytes);
}
HPX_ACTION(HPX_FUNCTION, 0,
           expansion_lco_init, expansion_lco_init_handler,
           HPX_POINTER, HPX_SIZE_T, HPX_POINTER, HPX_SIZE_T);


/// The set operation handler for the Expansion LCO
///
/// This takes one of three forms. The input to this is either a single integer
/// or a serialized Expansion. In the latter case, the reserved data at the
/// beginning of the expansion serialization is used to give the operation
/// code for the set.
void expansion_lco_operation_handler(ExpansionLCOHeader *lhs, void *rhs,
                                     size_t bytes) {
  int *code = static_cast<int *>(rhs);

  if (*code == kFinish) {
    assert(lhs->finished == 0);
    lhs->finished = 1;
  } else if (*code == kSchedule) {
    assert(lhs->finished == 0);
    lhs->scheduled += 1;
  } else if (*code == kContribute) {
    //create an expansion from the rhs
    char *input = static_cast<char *>(rhs);
    int *rhstype = reinterpret_cast<int *>(input + sizeof(int));
    int *n_digits = reinterpret_cast<int *>(input + sizeof(int) * 2);

    auto incoming = interpret_expansion(*rhstype, input, bytes,
                                        *n_digits);

    //create the expansion from the payload
    int *type = reinterpret_cast<int *>(lhs->payload + sizeof(int));
    auto local = interpret_expansion(*type, lhs->payload, lhs->payload_size,
                                     *n_digits);

    //add the one to the other
    local->add_expansion(incoming.get());

    //release the data, because these objects do not actually own those buffers
    local->release();
    incoming->release();

    //increment the counter
    lhs->arrived += 1;
    if (lhs->finished) {
      assert(lhs->arrived <= lhs->scheduled);
    }
  } else {
    assert(0 && "Incorrect code to expansion LCO");
  }
}
HPX_ACTION(HPX_FUNCTION, 0,
           expansion_lco_operation, expansion_lco_operation_handler,
           HPX_POINTER, HPX_POINTER, HPX_SIZE_T);


/// The predicate to detect triggering of the Expansion LCO
///
/// The expansion LCO is triggered if it has been finalized and the number
/// of contributions match the number of scheduled operations.
bool expansion_lco_predicate_handler(ExpansionLCOHeader *i, size_t bytes) {
  return (i->finished && (i->arrived == i->scheduled));
}
HPX_ACTION(HPX_FUNCTION, 0,
           expansion_lco_predicate, expansion_lco_predicate_handler,
           HPX_POINTER, HPX_SIZE_T);


/////////////////////////////////////////////////////////////////////
// HPX Stuff
/////////////////////////////////////////////////////////////////////


int expansion_s_to_m_handler(Source *sources, int n_src,
                             double cx, double cy, double cz, double scale,
                             hpx_addr_t expand, int type, int n_digits) {
  auto local = interpret_expansion(type, nullptr, 0, n_digits);
  auto multi =
    local->S_to_M(Point{cx, cy, cz}, sources, &sources[n_src], scale);
  size_t bytes = multi->bytes();
  char *serial = static_cast<char *>(multi->release());

  ExpansionRef total{type, expand, n_digits};
  total.contribute(bytes, serial);
  free(serial);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           expansion_s_to_m_action, expansion_s_to_m_handler,
           HPX_POINTER, HPX_INT, HPX_DOUBLE, HPX_DOUBLE, HPX_DOUBLE,
           HPX_DOUBLE, HPX_ADDR, HPX_INT, HPX_INT);


int expansion_s_to_l_handler(Source *sources, int n_src,
                             double cx, double cy, double cz, double scale,
                             hpx_addr_t expand, int type, int n_digits) {
  auto local = interpret_expansion(type, nullptr, 0, n_digits);
  auto multi = local->S_to_L(Point{cx, cy, cz}, sources, &sources[n_src], scale);
  size_t bytes = multi->bytes();
  char *serial = static_cast<char *>(multi->release());

  ExpansionRef total{type, expand, n_digits};
  total.contribute(bytes, serial);
  free(serial);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           expansion_s_to_l_action, expansion_s_to_l_handler,
           HPX_POINTER, HPX_INT, HPX_DOUBLE, HPX_DOUBLE, HPX_DOUBLE,
           HPX_DOUBLE, HPX_ADDR, HPX_INT, HPX_INT);


int expansion_m_to_m_handler(int type, hpx_addr_t expand, int n_digits,
                             int from_child, double s_size) {
  hpx_addr_t target = hpx_thread_current_target();
  //HACK: This action is local to the expansion, so we getref here with
  // whatever as the size and things are okay...
  ExpansionLCOHeader *ldata{nullptr};
  hpx_lco_getref(target, 1, (void **)&ldata);
  auto lexp = interpret_expansion(type, ldata->payload, ldata->payload_size,
                                  n_digits);
  auto translated = lexp->M_to_M(from_child, s_size);
  lexp->release();
  hpx_lco_release(target, ldata);

  size_t bytes = translated->bytes();
  void *transexpand = translated->release();

  ExpansionRef total{type, expand, n_digits};
  total.contribute(bytes, static_cast<char *>(transexpand));

  free(transexpand);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           expansion_m_to_m_action, expansion_m_to_m_handler,
           HPX_INT, HPX_ADDR, HPX_INT, HPX_INT, HPX_DOUBLE);


struct ExpansionMtoLParams {
  ExpansionRef total;
  Index s_index;
  double s_size;
  Index t_index;
  int n_digits;
};

int expansion_m_to_l_handler(ExpansionMtoLParams *parms, size_t UNUSED) {
  hpx_addr_t target = hpx_thread_current_target();
  //HACK: This action is local to the expansion, so we getref here with
  // whatever as the size and things are okay...
  ExpansionLCOHeader *ldata{nullptr};
  hpx_lco_getref(target, 1, (void **)&ldata);

  auto lexp = interpret_expansion(parms->total.type(), ldata->payload,
                                  ldata->payload_size, parms->n_digits);
  auto translated = lexp->M_to_L(parms->s_index, parms->s_size, parms->t_index);
  lexp->release();
  hpx_lco_release(target, ldata);

  size_t bytes = translated->bytes();
  void *transexpand = translated->release();

  parms->total.contribute(bytes, static_cast<char *>(transexpand));

  free(transexpand);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
           expansion_m_to_l_action, expansion_m_to_l_handler,
           HPX_POINTER, HPX_SIZE_T);


int expansion_l_to_l_handler(int type, hpx_addr_t expand, int n_digits,
                             int to_child, double t_size) {
  hpx_addr_t target = hpx_thread_current_target();
  //HACK: This action is local to the expansion, so we getref here with
  // whatever as the size and things are okay...
  ExpansionLCOHeader *ldata{nullptr};
  hpx_lco_getref(target, 1, (void **)&ldata);
  auto lexp = interpret_expansion(type, ldata->payload, ldata->payload_size,
                                  n_digits);
  auto translated = lexp->L_to_L(to_child, t_size);
  lexp->release();
  hpx_lco_release(target, ldata);

  size_t bytes = translated->bytes();
  void *transexpand = translated->release();

  ExpansionRef total{type, expand, n_digits};
  total.contribute(bytes, static_cast<char *>(transexpand));

  free(transexpand);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           expansion_l_to_l_action, expansion_l_to_l_handler,
           HPX_INT, HPX_ADDR, HPX_INT, HPX_INT, HPX_DOUBLE);


int expansion_m_to_t_handler(int n_targets, int n_digits, double scale,
                             hpx_addr_t targ, int type) {
  TargetRef targets{targ, n_targets};
  //HACK: This action is local to the expansion, so we getref here with
  // whatever as the size and things are okay...
  ExpansionLCOHeader *ldata{nullptr};
  hpx_lco_getref(hpx_thread_current_target(), 1, (void **)&ldata);
  targets.contribute_M_to_T(type, ldata->payload_size, ldata->payload,
                            n_digits, scale);
  hpx_lco_release(hpx_thread_current_target(), ldata);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           expansion_m_to_t_action, expansion_m_to_t_handler,
           HPX_INT, HPX_INT, HPX_DOUBLE, HPX_ADDR, HPX_INT);


int expansion_l_to_t_handler(int n_targets, int n_digits, double scale,
                             hpx_addr_t targ, int type) {
  TargetRef targets{targ, n_targets};
  //HACK: This action is local to the expansion, so we getref here with
  // whatever as the size and things are okay...
  ExpansionLCOHeader *ldata{nullptr};
  hpx_lco_getref(hpx_thread_current_target(), 1, (void **)&ldata);
  targets.contribute_L_to_T(type, ldata->payload_size, ldata->payload,
                            n_digits, scale);
  hpx_lco_release(hpx_thread_current_target(), ldata);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           expansion_l_to_t_action, expansion_l_to_t_handler,
           HPX_INT, HPX_INT, HPX_DOUBLE, HPX_ADDR, HPX_INT);


int expansion_s_to_t_handler(Source *sources, int n_sources, hpx_addr_t target,
                             int type) {
  TargetRef targets{target, 0};
  targets.contribute_S_to_T(type, n_sources, sources);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           expansion_s_to_t_action, expansion_s_to_t_handler,
           HPX_POINTER, HPX_INT, HPX_ADDR, HPX_INT);


int expansion_add_handler(hpx_addr_t expand, int type, int n_digits) {
  ExpansionRef total{type, expand, n_digits};

  ExpansionLCOHeader *ldata{nullptr};
  //HACK: This action is local to the expansion, so we getref here with
  // whatever as the size and things are okay...
  hpx_addr_t target = hpx_thread_current_target();
  hpx_lco_getref(target, 1, (void **)&ldata);
  total.contribute(ldata->payload_size, ldata->payload);
  hpx_lco_release(target, ldata);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           expansion_add_action, expansion_add_handler,
           HPX_ADDR, HPX_INT, HPX_INT);


int globalize_expansion_handler(void *payload, size_t bytes) {
  size_t total_size = sizeof(ExpansionLCOHeader) + bytes;
  hpx_addr_t gdata = hpx_lco_user_new(total_size, expansion_lco_init,
                                      expansion_lco_operation,
                                      expansion_lco_predicate, payload, bytes);
  assert(gdata != HPX_NULL);
  return HPX_THREAD_CONTINUE(gdata);
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
           globalize_expansion_action, globalize_expansion_handler,
           HPX_POINTER, HPX_SIZE_T);


/////////////////////////////////////////////////////////////////////
// Interface
/////////////////////////////////////////////////////////////////////


void ExpansionRef::destroy() {
  if (data_ != HPX_NULL) {
    hpx_lco_delete_sync(data_);
    data_ = HPX_NULL;
    type_ = 0;
  }
}


void ExpansionRef::S_to_M(Point center, SourceRef sources,
                          double scale) const {
  schedule();
  int nsrc = sources.n();
  double cx = center.x();
  double cy = center.y();
  double cz = center.z();
  hpx_call(sources.data(), expansion_s_to_m_action, HPX_NULL,
           &nsrc, &cx, &cy, &cz, &scale, &data_, &type_, &n_digits_);
}


void ExpansionRef::S_to_L(Point center, SourceRef sources,
                          double scale) const {
  schedule();
  int nsrc = sources.n();
  double cx = center.x();
  double cy = center.y();
  double cz = center.z();
  hpx_call(sources.data(), expansion_s_to_l_action, HPX_NULL,
           &nsrc, &cx, &cy, &cz, &scale, &data_, &type_, &n_digits_);
}


void ExpansionRef::M_to_M(ExpansionRef source, int from_child,
                          double s_size) const {
  assert(type_ == source.type());
  schedule();
  hpx_call_when(source.data(), source.data(), expansion_m_to_m_action, HPX_NULL,
                &type_, &data_, &n_digits_, &from_child, &s_size);
}


void ExpansionRef::M_to_L(ExpansionRef source, Index s_index, double s_size,
                          Index t_index) const {
  assert(type_ == source.type());
  schedule();
  ExpansionMtoLParams args{*this, s_index, s_size, t_index, n_digits_};
  hpx_call_when(source.data(), source.data(), expansion_m_to_l_action, HPX_NULL,
                &args, sizeof(args));
}


void ExpansionRef::L_to_L(ExpansionRef source, int to_child,
                          double t_size) const {
  assert(type_ == source.type());
  schedule();
  hpx_call_when(source.data(), source.data(), expansion_l_to_l_action, HPX_NULL,
                &type_, &data_, &n_digits_, &to_child, &t_size);
}


void ExpansionRef::M_to_T(TargetRef targets, double scale) const {
  targets.schedule(1);
  int nsend = targets.n();
  hpx_addr_t tsend = targets.data();
  hpx_call_when(data_, data_, expansion_m_to_t_action, HPX_NULL,
                &nsend, &n_digits_, &scale, &tsend, &type_);
}


void ExpansionRef::L_to_T(TargetRef targets, double scale) const {
  targets.schedule(1);
  int nsend = targets.n();
  hpx_addr_t tsend = targets.data();
  hpx_call_when(data_, data_, expansion_l_to_t_action, HPX_NULL,
                &nsend, &n_digits_, &scale, &tsend, &type_);
}


void ExpansionRef::S_to_T(SourceRef sources, TargetRef targets) const {
  targets.schedule(1);
  int n_src = sources.n();
  hpx_addr_t tsend = targets.data();
  hpx_call(sources.data(), expansion_s_to_t_action, HPX_NULL, &n_src, &tsend,
           &type_);
}


void ExpansionRef::add_expansion(ExpansionRef summand) {
  schedule();  //we are going to have another contribution
  hpx_call_when(summand.data(), summand.data(), expansion_add_action,
                HPX_NULL, &data_, &type_, &n_digits_);
}


void ExpansionRef::contribute(size_t bytes, char *payload) {
  int *code = reinterpret_cast<int *>(payload);
  *code = kContribute;
  hpx_lco_set_lsync(data_, bytes, payload, HPX_NULL);
}


std::unique_ptr<Expansion> ExpansionRef::get_new_expansion(Point center,
                                                           int n_digits) const {
  return create_expansion(type_, center, n_digits);
}


void ExpansionRef::finalize() const {
  if (data_ != HPX_NULL) {
    int code = kFinish;
    hpx_lco_set_lsync(data_, sizeof(code), &code, HPX_NULL);
  }
}


void ExpansionRef::schedule() const {
  if (data_ != HPX_NULL) {
    int code = kSchedule;
    hpx_lco_set_rsync(data_, sizeof(code), &code);
  }
}


//NOTE that this function takes ownership of the Expansion.
ExpansionRef globalize_expansion(std::unique_ptr<Expansion> exp,
                                 hpx_addr_t where) {
  if (exp == nullptr) {
    return ExpansionRef{0, HPX_NULL, -1};
  }

  //This is the init data for the LCO
  size_t bytes = exp->bytes();
  int type = exp->type();
  int n_digits = exp->accuracy();
  void *ldata = exp->release();

  hpx_addr_t retval{HPX_NULL};
  hpx_call_sync(where, globalize_expansion_action, &retval, sizeof(retval),
                ldata, bytes);
  free(ldata);

  return ExpansionRef{type, retval, n_digits};
}


} // namespace dashmm
