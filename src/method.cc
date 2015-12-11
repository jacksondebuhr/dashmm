/// \file src/method.cc
/// \brief Implementation of method

#include "include/method.h"

//C++

#include <hpx/hpx.h>

#include "include/methodref.h"


namespace dashmm {


/// The mapping from method type identifier into creation actions
///
/// This is local data that is replicated at each locality in the system.
std::map<int, hpx_action_t> *method_table_;


/// The first built-in method type identifier
constexpr int kFirstMethodType = 0;

/// The last built-in method type identifier
constexpr int kLastMethodType = 999;

constexpr int kFirstUserMethodType = 1000;
constexpr int kLastUserMethodType = 1999;


constexpr int kAllocateMethodTable = 0;
constexpr int kDeleteMethodTable = 1;


/////////////////////////////////////////////////////////////////////
// Actions
/////////////////////////////////////////////////////////////////////


int register_method_handler(int type, hpx_action_t creator, hpx_addr_t check) {
  int checkval = 0;
  if (method_table_->count(type) != 0) {
    int checkval = 1;
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
  }
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           manage_method_table_action, manage_method_table_handler,
           HPX_INT);


/////////////////////////////////////////////////////////////////////
// Internal routines
/////////////////////////////////////////////////////////////////////


std::unique_ptr<Method> create_method(int type, MethodSerial *data) {
  auto entry = method_table_->find(type);
  if (entry == method_table_->end()) {
    return std::unique_ptr<Method>{nullptr};
  }
  method_creation_function_t func =
      reinterpret_cast<method_creation_function_t>(
        hpx_action_get_handler(entry->second)
      );
  return std::unique_ptr<Method>{func(sizeof(MethodSerial) + data->size, data)};
}


void init_method_table() {
  int input = kAllocateMethodTable;
  hpx_bcast_rsync(manage_method_table_action, &input);
}


void fini_method_table() {
  int input = kDeleteMethodTable;
  hpx_bcast_rsync(manage_method_table_action, &input);
}


/////////////////////////////////////////////////////////////////////
// In the public interface to dashmm
/////////////////////////////////////////////////////////////////////


ReturnCode register_method(int type, hpx_action_t creator) {
  if (type < kFirstUserMethodType || type > kLastUserMethodType) {
    return kDomainError;
  }

  int nlocs = hpx_get_num_ranks();
  hpx_addr_t checker = hpx_lco_reduce_new(nlocs, sizeof(int),
                                          int_sum_ident_op, int_sum_op);
  assert(checker != HPX_NULL);
  hpx_bcast_lsync(register_method_action, HPX_NULL, &type, &creator, &checker);
  int checkval{0};
  hpx_lco_get(checker, sizeof(int), &checkval);
  return (checkval == 0 ? kSuccess : kDomainError);
}


} // namespace dashmm
