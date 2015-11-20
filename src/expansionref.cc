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
    char *payload = static_cast<char *>(i) + sizeof(ExpansionLCOHeader);
    auto local = interpret_expansion(payload, head->payload_size);

    //create an expansion from the rhs
    char *input = static_cast<char *>(rhs) + sizeof(int);
    auto incoming = interpret_expansion(input, bytes - sizeof(int));

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
  setup_local_expansion();
  return exp_->S_to_M(center, first, last);
}


//TODO
std::unique_ptr<Expansion> ExpansionRef::S_to_L(Point center,
                                      Source *first, Source *last) const {
  setup_local_expansion();
  return exp_->S_to_L(center, first, last);
}


//TODO
std::unique_ptr<Expansion> ExpansionRef::M_to_M(int from_child,
                                                double s_size) const {
  setup_local_expansion();
  return exp_->M_to_M(from_child, s_size);
}


//TODO
std::unique_ptr<Expansion> ExpansionRef::M_to_L(Index s_index, double s_size,
                                  Index t_index) const {
  setup_local_expansion();
  return exp_->M_to_L(s_index, s_size, t_index);
}


//TODO
std::unique_ptr<Expansion> ExpansionRef::L_to_L(int to_child,
                                                double t_size) const {
  setup_local_expansion();
  return exp_->L_to_L(to_child, t_size);
}


//TODO
void ExpansionRef::M_to_T(Target *first, Target *last) const {
  setup_local_expansion();
  exp_->M_to_T(first, last);
}


//TODO
void ExpansionRef::L_to_T(Target *first, Target *last) const {
  setup_local_expansion();
  exp_->L_to_T(first, last);
}


//TODO
void ExpansionRef::S_to_T(Source *s_first, Source *s_last,
                          Target *t_first, Target *t_last) const {
  setup_local_expansion();
  exp_->S_to_T(s_first, s_last, t_first, t_last);
}


//TODO
void ExpansionRef::add_expansion(const Expansion *temp1) {
  setup_local_expansion();
  exp_->add_expansion(temp1);
  save_to_global();
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
  void *data = exp->release();
  int *type = static_cast<int *>(data);

  size_t total_size = sizeof(ExpansionLCOHeader) + bytes;

  hpx_addr_t data = hpx_lco_user_new(total_size, expansion_lco_init,
                                     expansion_lco_operation,
                                     expansion_lco_predicate, data, bytes);
  assert(data != HPX_NULL);

  return ExpansionRef{*type, data};
}


} // namespace dashmm
