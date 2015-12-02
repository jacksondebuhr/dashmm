#include "include/expansionref.h"

#include <cassert>

#include <hpx/hpx.h>

//other dashmm


namespace dashmm {


/////////////////////////////////////////////////////////////////////
// Stuff for the user LCO
/////////////////////////////////////////////////////////////////////


struct ExpansionLCOHeader {
  int arrived;
  int scheduled;
  int finished;
  size_t payload_size;
};

enum ExpansionLCOSetCodes {
  kFinish = 1,
  kSchedule = 2,
  kContribute = 3
};


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


void expansion_lco_operation_handler(void *lhs, void *rhs, size_t bytes) {
  int *code = static_cast<int *>(rhs);
  if (*code == kFinish) {
    ExpansionLCOHeader *head = static_cast<ExpansionLCOHeader *>(lhs);
    head->finished = 1;
  } else if (*code == kSchedule) {
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


bool expansion_lco_predicate_handler(ExpansionLCOHeader *i, size_t bytes) {
  return (i->finished && (i->arrived == i->scheduled));
}
HPX_ACTION(HPX_FUNCTION, 0,
           expansion_lco_predicate, expansion_lco_predicate_handler);


/////////////////////////////////////////////////////////////////////
// HPX Stuff
/////////////////////////////////////////////////////////////////////


int expansion_m_to_t_handler(int n_targets, hpx_addr_t targ, int type) {
  TargetRef targets{targ, n_targets};
  //HACK: This action is local to the expansion, so we getref here with
  // whatever as the size and things are okay...
  ExpansionLCOHeader *ldata{nullptr};
  hpx_lco_getref(hpx_thread_current_target(), 1, &ldata);
  char *payload = static_cast<char *>(ldata) + sizeof(ExpansionLCOHeader);
  targets.contribute_M_to_T(type, ldata->payload_size, payload);

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


//TODO
std::unique_ptr<Expansion> ExpansionRef::S_to_M(Point center,
                                      Source *first, Source *last) const {
  //
}


//TODO
std::unique_ptr<Expansion> ExpansionRef::S_to_L(Point center,
                                      Source *first, Source *last) const {
  //
}


//TODO
std::unique_ptr<Expansion> ExpansionRef::M_to_M(int from_child,
                                                double s_size) const {
  //
}


//TODO
std::unique_ptr<Expansion> ExpansionRef::M_to_L(Index s_index, double s_size,
                                  Index t_index) const {
  //
}


//TODO
std::unique_ptr<Expansion> ExpansionRef::L_to_L(int to_child,
                                                double t_size) const {
  //
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


void ExpansionRef::add_expansion(std::unique_ptr<Expansion> summand) {
  schedule();  //we are going to have another contribution
  size_t bytes = summand->bytes();
  char *payload = summand->release();
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
