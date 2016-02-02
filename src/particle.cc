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


/// \file src/particle.cc
/// \brief Source and Target object implementations


#include "include/particle.h"

#include <cassert>
#include <cstring>

#include <hpx/hpx.h>

#include "include/expansion.h"


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
  int code;
  int type;
  int count;
  Source sources[];
};

struct TargetRefLCOSetMtoTData {
  int code;
  int type;
  int n_digits;
  double scale;
  size_t bytes;
  char data[];
};

struct TargetRefLCOSetLtoTData {
  int code;
  int type;
  int n_digits;
  double scale;
  size_t bytes;
  char data[];
};

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
           targetref_lco_init, targetref_lco_init_handler,
           HPX_POINTER, HPX_SIZE_T, HPX_POINTER, HPX_SIZE_T);


void targetref_lco_operation_handler(TargetRefLCOData *lhs,
                                     void *rhs, size_t bytes) {
  int *code = static_cast<int *>(rhs);
  if (*code == kSetOnly) {    //this is a pair of ints, a code and a count
    lhs->scheduled += code[1];
  } else if (*code == kStoT) {
    TargetRefLCOSetStoTData *input =
        static_cast<TargetRefLCOSetStoTData *>(rhs);
    auto expansion = interpret_expansion(input->type, nullptr, 0, -1);
    expansion->S_to_T(input->sources, &input->sources[input->count],
                      lhs->targets, &lhs->targets[lhs->count]);
    expansion->release();
    lhs->arrived += 1;
  } else if (*code == kMtoT) {
    TargetRefLCOSetMtoTData *input =
        static_cast<TargetRefLCOSetMtoTData *>(rhs);
    auto expansion = interpret_expansion(input->type, input->data,
                                         input->bytes, input->n_digits);
    expansion->M_to_T(lhs->targets, &lhs->targets[lhs->count], input->scale);
    expansion->release();
    lhs->arrived += 1;
  } else if (*code == kLtoT) {
    TargetRefLCOSetLtoTData *input =
        static_cast<TargetRefLCOSetLtoTData *>(rhs);
    auto expansion = interpret_expansion(input->type, input->data,
                                         input->bytes, input->n_digits);
    expansion->L_to_T(lhs->targets, &lhs->targets[lhs->count], input->scale);
    expansion->release();
    lhs->arrived += 1;
  } else if (*code == kFinish) {
    lhs->finished = 1;
  } else {
    assert(0 && "Incorrect code to TargetRef LCO");
  }
}
HPX_ACTION(HPX_FUNCTION, 0,
           targetref_lco_operation, targetref_lco_operation_handler,
           HPX_POINTER, HPX_POINTER, HPX_SIZE_T);


bool targetref_lco_predicate_handler(TargetRefLCOData *i, size_t bytes) {
  return i->finished && (i->arrived == i->scheduled);
}
HPX_ACTION(HPX_FUNCTION, 0,
           targetref_lco_predicate, targetref_lco_predicate_handler,
           HPX_POINTER, HPX_SIZE_T);


/////////////////////////////////////////////////////////////////////
// Interface
/////////////////////////////////////////////////////////////////////


void SourceRef::destroy() {
  if (data_ != HPX_NULL) {
    hpx_gas_free_sync(data_);
    data_ = HPX_NULL;
    n_ = 0;
  }
}


SourceRef SourceRef::slice(size_t offset, size_t n) const {
  if (offset > n_) {
    return SourceRef{HPX_NULL, 0, 0, 0, 0, 0};
  }
  if (offset + n > n_) {
    return SourceRef{HPX_NULL, 0, 0, 0, 0, 0};
  }
  hpx_addr_t addr = hpx_addr_add(data_, sizeof(Source) * offset,
                                 sizeof(Source) * n_tot_);
  return SourceRef{addr, n, n_tot_, record_size_, pos_offset_, q_offset_};
}


TargetRef::TargetRef(Target *targets, int n) {
  size_t init_size = sizeof(TargetRefLCOInitData) + sizeof(Target) * n;
  TargetRefLCOInitData *init =
      reinterpret_cast<TargetRefLCOInitData *>(new char [init_size]);
  assert(init);
  init->count = n;
  memcpy(init->targets, targets, sizeof(Target) * n);

  size_t total_size = sizeof(TargetRefLCOData) + sizeof(Target) * n;

  data_ = hpx_lco_user_new(total_size, targetref_lco_init,
                           targetref_lco_operation, targetref_lco_predicate,
                           init, init_size);
  assert(data_ != HPX_NULL);
  delete [] init;
  n_ = n;
}


void TargetRef::destroy() {
  if (data_ != HPX_NULL) {
    hpx_lco_delete_sync(data_);
    data_ = HPX_NULL;
    n_ = 0;
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


void TargetRef::contribute_S_to_T(int type, int n, Source *sources) const {
  //NOTE: We assume this is called local to the sources
  size_t inputsize = sizeof(TargetRefLCOSetStoTData)
                     + sizeof(Source) * n;
  TargetRefLCOSetStoTData *input =
      reinterpret_cast<TargetRefLCOSetStoTData *>(new char [inputsize]);
  assert(input);
  input->code = kStoT;
  input->type = type;
  input->count = n;
  memcpy(input->sources, sources, sizeof(Source) * n);

  //call set with appropriate data
  hpx_lco_set_lsync(data_, inputsize, input, HPX_NULL);

  delete [] input;
}


void TargetRef::contribute_M_to_T(int type, size_t bytes, void *data,
                                  int n_digits, double scale) const {
  size_t inputsize = sizeof(TargetRefLCOSetMtoTData) + bytes;
  TargetRefLCOSetMtoTData *input =
      reinterpret_cast<TargetRefLCOSetMtoTData *>(new char [inputsize]);
  assert(input);
  input->code = kMtoT;
  input->type = type;
  input->n_digits = n_digits;
  input->scale = scale;
  input->bytes = bytes;
  memcpy(input->data, data, bytes);
  hpx_lco_set_lsync(data_, inputsize, input, HPX_NULL);
  delete [] input;
}


void TargetRef::contribute_L_to_T(int type, size_t bytes, void *data,
                                  int n_digits, double scale) const {
  size_t inputsize = sizeof(TargetRefLCOSetLtoTData) + bytes;
  TargetRefLCOSetLtoTData *input =
      reinterpret_cast<TargetRefLCOSetLtoTData *>(new char [inputsize]);
  assert(input);
  input->code = kLtoT;
  input->type = type;
  input->n_digits = n_digits;
  input->scale = scale;
  input->bytes = bytes;
  memcpy(input->data, data, bytes);
  hpx_lco_set_lsync(data_, inputsize, input, HPX_NULL);
}


} // namespace dashmm
