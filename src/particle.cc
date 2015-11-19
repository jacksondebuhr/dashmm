#include "include/particle.h"

// C/C++

#include <hpx/hpx.h>

// other dashmm


namespace dashmm {


/////////////////////////////////////////////////////////////////////
// Stuff for the user LCO
/////////////////////////////////////////////////////////////////////


struct TargetRefLCOData {
  int arrived;
  int scheduled;
  int finished;
  int count;
  Target targets[];
};

struct TargetRefLCOInitData {
  int count;
  Target targets[];
};

struct TargetRefLCOSetStoTData {
  int type;
  int count;
  Source sources[];
}

struct TargetRefLCOSetMtoTData {
  int type;
  size_t bytes;
  char data[];
}

struct TargetRefLCOSetLtoTData {
  int type;
  size_t bytes;
  char data[];
}

enum TargetRefLCOSetCodes {
  kSetOnly = 1,
  kStoT = 2,
  kMtoT = 3,
  kLtoT = 4,
  kFinish = 5
};


void targetref_lco_init_handler(TargetRefLCOData *i, size_t bytes,
                                TargetRefLCOInitData *init, size_t init_bytes) {
  i->arrived = 0;
  i->scheduled = 0;
  i->finished = 0;
  i->count = init->count;
  memcpy(i->targets, init->targets, sizeof(Target) * init->count);
}
HPX_ACTION(HPX_FUNCTION, 0,
           targetref_lco_init, targetref_lco_init_handler);


void targetref_lco_operation_handler(TargetRefLCOData *lhs,
                                     void *rhs, size_t bytes) {
  int *code = static_cast<int *>(rhs);
  if (*code == kSetOnly) {    //this is a pair of ints, a code and a count
    lhs->arrived += code[1];
  } else if (*code == kStoT) {
    TargetRefLCOSetStoTData *input =
        static_cast<TargetRefLCOSetStoTData *>(rhs);
    auto expansion = interpret_expansion(input->type, nullptr, 0);
    expansion->S_to_T(input->sources, &input->sources[input->count],
                      lhs->targets, &lhs->targets[lhs->count]);
    lhs->arrived += 1;
  } else if (*code == kMtoT) {
    TargetRefLCOSetMtoTData *input =
        static_cast<TargetRefLCOSetMtoTData *>(rhs);
    auto expansion = interpret_expansion(input->type, input->data,
                                         input->bytes);
    expansion->M_to_T(lhs->targets, &lhs->targets[lhs->count]);
    lhs->arrived += 1;
  } else if (*code == kLtoT) {
    TargetRefLCOSetLtoTData *input =
        static_cast<TargetRefLCOSetLtoTData *>(rhs);
    auto expansion = interpret_expansion(input->type, input->data,
                                         input->bytes);
    expansion->L_to_T(lhs->targets, &lhs->targets[lhs->count]);
    lhs->arrived += 1;
  } else if (*code == kFinish) {
    lhs->finished = 1;
  } else {
    assert(0 && "Incorrect code to TargetRef LCO");
  }
}
HPX_ACTION(HPX_FUNCTION, 0,
           targetref_lco_operation, targetref_lco_operation_handler);


bool targetref_lco_predicate_handler(TargetRefLCOData *i, size_t bytes) {
  return i->finished && (i->arrived == i->scheduled);
}
HPX_ACTION(HPX_FUNCTION, 0,
           targetref_lco_predicate, targetref_lco_predicate_handler);


/////////////////////////////////////////////////////////////////////
// Interface
/////////////////////////////////////////////////////////////////////


SourceRef::SourceRef(Source *sources, int n) {
  data_ = hpx_gas_alloc_local(1, sizeof(Source) * n, 0);
  assert(data_ != HPX_NULL);
  n_ = n;
  Source *local{nullptr};
  assert(hpx_gas_try_pin(data_, (void **)local));
  memcpy(local, sources, sizeof(Source) * n);
  hpx_gas_unpin(data_);
}


void SourceRef::destroy() {
  if (data_ != HPX_NULL) {
    hpx_gas_free_sync(data_);
    data_ = HPX_NULL;
    n = 0;
  }
}


TargetRef::TargetRef(Target *targets, int n) {
  size_t init_size = sieof(TargetRefLCOInitData) + sizeof(Target) * n;
  TargetRefLCOInitData *init = malloc(init_size);
  assert(init);
  init->count = n;
  memcpy(init->targets, targets, sizeof(Target) * n);

  data_ = hpx_lco_user_new(init_size, targetref_lco_init,
                           targetref_lco_operation, targetref_lco_predicate,
                           init, init_size);
  assert(data_ != HPX_NULL);
  n_ = n;
}


void TargetRef::destroy() {
  if (data_ != HPX_NULL) {
    hpx_lco_delete_sync(data_);
    data_ = HPX_NULL;
    n = 0;
  }
}


void TargetRef::finalize() const {
  if (data_ != HPX_NULL) {
    int code = kFinish;
    hpx_lco_set_lsync(data_, sizeof(code), &code, HPX_NULL);
  }
}


//NOTE: This one must be remote complete so that there is no timing issue
// between scheduling an input and finalizing the schedule.
void TargetRef::schedule(int num) const {
  if (data_ != HPX_NULL) {
    int input[2]{kSetOnly, num};
    hpx_lco_set_rsync(data_, sizeof(int) * 2, input);
  }
}


} // namespace dashmm