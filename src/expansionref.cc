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
void expansion_lco_init_handler(void *i, size_t bytes,
                                void *init, size_t init_bytes) {
  ExpansionLCOHeader *head = static_cast<ExpansionLCOHeader *>(i);
  i->arrived = 0;
  i->scheduled = 0;
  i->finished = 0;
  i->payload_size = init_bytes;
  char *payload = static_cast<char *>(i) + sizeof(ExpansionLCOHeader);
  memcpy(payload, init, init_bytes);
}
HPX_ACTION(HPX_FUNCTION, 0,
           expansion_lco_init, expansion_lco_init_handler);


/// The set operation handler for the Expansion LCO
///
/// This takes one of three forms. The input to this is either a single integer
/// or a serialized Expansion. In the latter case, the reserved data at the
/// beginning of the expansion serialization is used to give the operation
/// code for the set.
void expansion_lco_operation_handler(void *lhs, void *rhs, size_t bytes) {
  int *code = static_cast<int *>(rhs);
  if (*code == kFinish) {
    ExpansionLCOHeader *head = static_cast<ExpansionLCOHeader *>(lhs);
    head->finished = 1;
  } else if (*code == kSchedule) {
    assert(head->finished == 0);
    ExpansionLCOHeader *head = static_cast<ExpansionLCOHeader *>(lhs);
    head->scheduled += 1;
  } else if (*code == kContribute) {
    //create the expansion from the payload
    ExpansionLCOHeader *head = static_cast<ExpansionLCOHeader *>(lhs);
    char *payload = static_cast<char *>(lhs) + sizeof(ExpansionLCOHeader);
    int *type = static_cast<int *>(payload + sizeof(int));
    auto local = interpret_expansion(*type, payload, head->payload_size);

    //create an expansion from the rhs
    char *input = static_cast<char *>(rhs);
    int *rhstype = static_cast<int *>(input + sizeof(int));
    auto incoming = interpret_expansion(*rhstype, input, bytes - sizeof(int));

    //add the one to the other
    local->add_expansion(incoming.get());

    //release the data, because these objects do not actually own those buffers
    local->release();
    incoming->release();

    //increment the counter
    head->arrived += 1;
  } else {
    assert(0 && "Incorrect code to expansion LCO");
  }
}
HPX_ACTION(HPX_FUNCTION, 0,
           expansion_lco_operation, expansion_lco_operation_handler);


/// The predicate to detect triggering of the Expansion LCO
///
/// The expansion LCO is triggered if it has been finalized() and the number
/// of contributions match the number of scheduled operations.
bool expansion_lco_predicate_handler(ExpansionLCOHeader *i, size_t bytes) {
  return (i->finished && (i->arrived == i->scheduled));
}
HPX_ACTION(HPX_FUNCTION, 0,
           expansion_lco_predicate, expansion_lco_predicate_handler);


/////////////////////////////////////////////////////////////////////
// HPX Stuff
/////////////////////////////////////////////////////////////////////


int expansion_s_to_m_handler(Source *sources, int n_src, double cx, double cy,
                             double cz, hpx_addr_t expand, int type) {
  auto local = interpret_expansion(type, nullptr, 0);
  auto multi = local->S_to_M(Point{cx, cy, cz}, sources, &sources[n_src]);
  size_t bytes = multi->bytes();
  char *serial = static_cast<char *>(multi->release());

  ExpansionRef total{type, expand};
  total.contribute(bytes, serial);
  free(serial);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           expansion_s_to_m_action, expansion_s_to_m_handler,
           HPX_POINTER, HPX_INT, HPX_DOUBLE, HPX_DOUBLE, HPX_DOUBLE,
           HPX_ADDR_T, HPX_INT);


int expansion_s_to_l_handler(Source *sources, int n_src, double cx, double cy,
                             double cz, hpx_addr_t expand, int type) {
  auto local = interpret_expansion(type, nullptr, 0);
  auto multi = local->S_to_L(Point{cx, cy, cz}, sources, &sources[n_src]);
  size_t bytes = multi->bytes();
  char *serial = static_cast<char *>(multi->release());

  ExpansionRef total{type, expand};
  total.contribute(bytes, serial);
  free(serial);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           expansion_s_to_m_action, expansion_s_to_m_handler,
           HPX_POINTER, HPX_INT, HPX_DOUBLE, HPX_DOUBLE, HPX_DOUBLE,
           HPX_ADDR_T, HPX_INT);


int expansion_m_to_m_handler(int type, hpx_addr_t expand, int from_child,
                             double s_size) {
  hpx_addr_t target = hpx_thread_current_target();
  //HACK: This action is local to the expansion, so we getref here with
  // whatever as the size and things are okay...
  ExpansionLCOHeader *ldata{nullptr};
  hpx_lco_getref(target, 1, &ldata);
  char *payload = static_cast<char *>(ldata) + sizeof(ExpansionLCOHeader);
  auto lexp = interpret_expansion(type, payload, ldata->payload_size);
  auto translated = lexp->M_to_M(from_child, s_size);
  lexp->release();
  hpx_lco_release(target, ldata);

  size_t bytes = translated->bytes();
  void *transexpand = translated->release();

  ExpansionRef total{type, expand};
  total.contribute(bytes, static_cast<char *>(transexpand));

  free(transexpand);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           expansion_m_to_m_action, expansion_m_to_m_handler,
           HPX_INT, HPX_ADDR_T, HPX_INT, HPX_DOUBLE);


struct ExpansionMtoLParams {
  ExpansionRef total;
  Index s_index;
  double s_size;
  Index t_index;
};

int expansion_m_to_l_action(ExpansionMtoLParams *parms, size_t UNUSED) {
  hpx_addr_t target = hpx_thread_current_target();
  //HACK: This action is local to the expansion, so we getref here with
  // whatever as the size and things are okay...
  ExpansionLCOHeader *ldata{nullptr};
  hpx_lco_getref(target, 1, &ldata);
  char *payload = static_cast<char *>(ldata) + sizeof(ExpansionLCOHeader);
  auto lexp = interpret_expansion(parms->total.type(), payload,
                                  ldata->payload_size);
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


int expansion_l_to_l_handler(int type, hpx_addr_t expand, int to_child,
                             double t_size) {
  hpx_addr_t target = hpx_thread_current_target();
  //HACK: This action is local to the expansion, so we getref here with
  // whatever as the size and things are okay...
  ExpansionLCOHeader *ldata{nullptr};
  hpx_lco_getref(target, 1, &ldata);
  char *payload = static_cast<char *>(ldata) + sizeof(ExpansionLCOHeader);
  auto lexp = interpret_expansion(type, payload, ldata->payload_size);
  auto translated = lexp->L_to_L(to_child, t_size);
  lexp->release();
  hpx_lco_release(target, ldata);

  size_t bytes = translated->bytes();
  void *transexpand = translated->release();

  ExpansionRef total{type, expand};
  total.contribute(bytes, static_cast<char *>(transexpand));

  free(transexpand);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           expansion_l_to_l_action, expansion_l_to_l_handler,
           HPX_INT, HPX_ADDR_T, HPX_INT, HPX_DOUBLE);


int expansion_m_to_t_handler(int n_targets, hpx_addr_t targ, int type) {
  TargetRef targets{targ, n_targets};
  //HACK: This action is local to the expansion, so we getref here with
  // whatever as the size and things are okay...
  ExpansionLCOHeader *ldata{nullptr};
  hpx_lco_getref(hpx_thread_current_target(), 1, &ldata);
  char *payload = static_cast<char *>(ldata) + sizeof(ExpansionLCOHeader);
  targets.contribute_M_to_T(type, ldata->payload_size, payload);
  hpx_lco_release(hpx_thread_current_target(), ldata);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           expansion_m_to_t_action, expansion_m_to_t_handler,
           HPX_INT, HPX_ADDR_T, HPX_INT);


int expansion_l_to_t_handler(int n_targets, hpx_addr_t targ, int type) {
  TargetRef targets{targ, n_targets};
  //HACK: This action is local to the expansion, so we getref here with
  // whatever as the size and things are okay...
  ExpansionLCOHeader *ldata{nullptr};
  hpx_lco_getref(hpx_thread_current_target(), 1, &ldata);
  char *payload = static_cast<char *>(ldata) + sizeof(ExpansionLCOHeader);
  targets.contribute_L_to_T(type, ldata->payload_size, payload);
  hpx_lco_release(hpx_thread_current_target(), ldata);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           expansion_l_to_t_action, expansion_l_to_t_handler,
           HPX_INT, HPX_ADDR_T, HPX_INT);


int expansion_s_to_t_handler(Source *sources, int n_sources, hpx_addr_t target,
                             int type) {
  TargetRef targets{target, 0};
  targets.contribute_S_to_T(type, n_sources, sources);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           expansion_s_to_t_action, expansion_s_to_t_handler,
           HPX_POINTER, HPX_INT, HPX_ADDR_T, HPX_INT);


int expansion_add_handler(hpx_addr_t expand, int type) {
  ExpansionRef total{type, expand};

  ExpansionLCOHeader *ldata{nullptr};
  //HACK: This action is local to the expansion, so we getref here with
  // whatever as the size and things are okay...
  hpx_addr_t target = hpx_thread_current_target();
  hpx_lco_getref(target, 1, &ldata);
  char *payload = static_cast<char *>(ldata) + sizeof(ExpansionLCOHeader);
  total.contribute(ldata->payload_size, payload);
  hpx_lco_release(target, ldata);

  hpx_lco_delete_sync(target);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           expansion_add_action, expansion_add_handler,
           HPX_ADDR_T, HPX_INT);


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


void ExpansionRef::S_to_M(Point center, SourceRef sources) const {
  schedule();
  int nsrc = sources.n();
  double cx = center.x();
  double cy = center.y();
  double cz = center.z();
  hpx_call(sources.data(), expansion_s_to_m_action, HPX_NULL,
           &nsrc, &cx, &cy, &cz, &data_, &type_);
}


void ExpansionRef::S_to_L(Point center, SourceRef sources) const {
  schedule();
  int nsrc = sources.n();
  double cx = center.x();
  double cy = center.y();
  double cz = center.z();
  hpx_call(sources.data(), expansion_s_to_l_action, HPX_NULL,
           &nsrc, &cx, &cy, &cz, &data_, &type_);
}


void ExpansionRef::M_to_M(ExpansionRef source, int from_child,
                          double s_size) const {
  assert(type_ == source.type());
  schedule();
  hpx_call_when(source.data(), source.data(), expansion_m_to_m_action, HPX_NULL,
                &type_, &data_, &from_child, &s_size);
}


void ExpansionRef::M_to_L(ExpansionRef source, Index s_index, double s_size,
                          Index t_index) const {
  assert(type_ == source.type());
  schedule();
  ExpansionMtoLParams args{*this, s_index, s_size, t_index};
  hpx_call_when(source.data(), source.data(), expansion_m_to_l_action, HPX_NULL,
                &args, sizeof(args));
}


void ExpansionRef::L_to_L(ExpansionRef source, int to_child,
                          double t_size) const {
  assert(type_ == source.type());
  schedule();
  hpx_call_when(source.data(), source.data(), expansion_l_to_l_action, HPX_NULL,
                &type_, &data_, &to_child, &t_size);
}


void ExpansionRef::M_to_T(TargetRef targets) const {
  targets.schedule(1);
  int nsend = targets.n();
  hpx_addr_t tsend = targets.data();
  hpx_call_when(data_, data_, expansion_m_to_t_action, HPX_NULL,
                &nsend, &tsend, &type_);
}


void ExpansionRef::L_to_T(TargetRef targets) const {
  targets.schedule(1);
  int nsend = targets.n();
  hpx_addr_t tsend = targets.data();
  hpx_call_when(data_, data_, expansion_l_to_t_action, HPX_NULL,
                &nsend, &tsend, &type_);
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
                HPX_NULL, &data_, &type_);
}


void ExpansionRef::contribute(size_t bytes, char *payload) {
  int *code = static_cast<int *>(payload);
  *code = kContribute;
  hpx_lco_set_lsync(data_, bytes, payload, HPX_NULL);
}


std::unique_ptr<Expansion> ExpansionRef::get_new_expansion(Point center) const {
  return create_expansion(type_, center);
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
ExpansionRef globalize_expansion(std::unique_ptr<Expansion> exp) {
  if (exp == nullptr) {
    return ExpansionRef{0, HPX_NULL};
  }

  //This is the init data for the LCO
  size_t bytes = exp->bytes();
  void *ldata = exp->release();
  char *offset = static_cast<char *>(ldata);
  int *ptype = static_cast<int *>(offset + sizeof(int));
  int type = *ptype;

  size_t total_size = sizeof(ExpansionLCOHeader) + bytes;

  hpx_addr_t gdata = hpx_lco_user_new(total_size, expansion_lco_init,
                                      expansion_lco_operation,
                                      expansion_lco_predicate, data, bytes);
  assert(gdata != HPX_NULL);
  free(ldata);

  return ExpansionRef{type, data};
}


} // namespace dashmm
