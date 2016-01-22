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


/// \file src/method.cc
/// \brief Implementation of method

#include "include/method.h"

#include <map>

#include <hpx/hpx.h>

#include "include/ids.h"
#include "include/methodref.h"
#include "include/reductionops.h"


namespace dashmm {


/// The mapping from method type identifier into creation actions
///
/// This is local data that is replicated at each locality in the system.
std::map<int, hpx_action_t> *method_table_;


constexpr int kAllocateMethodTable = 0;
constexpr int kDeleteMethodTable = 1;


/////////////////////////////////////////////////////////////////////
// Actions
/////////////////////////////////////////////////////////////////////


int register_method_handler(int type, hpx_action_t creator, hpx_addr_t check) {
  int checkval = 0;
  if (method_table_->count(type) != 0) {
    checkval = 1;
  } else {
    (*method_table_)[type] = creator;
  }
  hpx_lco_set_lsync(check, sizeof(int), &checkval, HPX_NULL);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           register_method_action, register_method_handler,
           HPX_INT, HPX_ACTION_T, HPX_ADDR);


int manage_method_table_handler(int op) {
  if (op == kAllocateMethodTable) {
    method_table_ = new std::map<int, hpx_action_t>{};
  } else if (op == kDeleteMethodTable) {
    delete method_table_;
    method_table_ = nullptr;
  }
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           manage_method_table_action, manage_method_table_handler,
           HPX_INT);


/////////////////////////////////////////////////////////////////////
// Internal routines
/////////////////////////////////////////////////////////////////////


ReturnCode register_method(int type, hpx_action_t creator, int user) {
  if (user) {
    if (type < kFirstUserMethodType || type > kLastUserMethodType) {
      return kDomainError;
    }
  } else {
    if (type < kFirstMethodType || type > kLastMethodType) {
      return kDomainError;
    }
  }

  int nlocs = hpx_get_num_ranks();
  hpx_addr_t checker = hpx_lco_reduce_new(nlocs, sizeof(int),
                                          int_sum_ident_op, int_sum_op);
  assert(checker != HPX_NULL);
  hpx_bcast_lsync(register_method_action, HPX_NULL, &type, &creator, &checker);
  int checkval{0};
  hpx_lco_get(checker, sizeof(int), &checkval);
  hpx_lco_delete_sync(checker);
  return (checkval == 0 ? kSuccess : kDomainError);
}


void init_method_table() {
  int input = kAllocateMethodTable;
  hpx_bcast_rsync(manage_method_table_action, &input);
}


void fini_method_table() {
  int input = kDeleteMethodTable;
  hpx_bcast_rsync(manage_method_table_action, &input);
}


Method *create_method(int type, MethodSerial *data) {
  auto entry = method_table_->find(type);
  if (entry == method_table_->end()) {
    return nullptr;
  }
  method_creation_function_t func =
      reinterpret_cast<method_creation_function_t>(
        hpx_action_get_handler(entry->second)
      );
  return func(sizeof(MethodSerial) + data->size, data);
}


/////////////////////////////////////////////////////////////////////
// In the public interface to dashmm
/////////////////////////////////////////////////////////////////////


ReturnCode register_user_method(int type, hpx_action_t creator) {
  return register_method(type, creator, 1);
}


} // namespace dashmm
