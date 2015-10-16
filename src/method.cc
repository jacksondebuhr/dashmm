#include "include/method.h"

//C++

#include <hpx/hpx.h>

#include "include/methodref.h"


namespace dashmm {


//The mapping from type to table entries
std::map<int, hpx_action_t> method_table_;


constexpr int kFirstMethodType = 0;
constexpr int kLastMethodType = 999;
constexpr int kFirstUserMethodType = 1000;
constexpr int kLastUserMethodType = 1999;


/////////////////////////////////////////////////////////////////////
// Actions
/////////////////////////////////////////////////////////////////////


int register_method_handler(int type, hpx_action_t creator, hpx_addr_t check) {
  int checkval = 0;
  if (method_table_[type].count() != 0) {
    int checkval = 1;
  } else {
    method_table_[type] = creator;
  }
  hpx_lco_set_lsync(check, sizeof(int), &checkval, HPX_NULL);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           register_method_action, register_method_handler,
           HPX_INT, HPX_ACTION_T, HPX_ADDR_T);


/////////////////////////////////////////////////////////////////////
// Internal routines
/////////////////////////////////////////////////////////////////////


std::unique_ptr<Method> create_method(int type, size_t size, void *data) {
  auto entry = method_table_.find(type);
  if (entry == method_table_.end()) {
    return std::unique_ptr<Method>{nullptr};
  }
  method_creation_function_t func =
      reinterpret_cast<method_creation_function_t>(
        hpx_action_get_handler(entry->second)
      );
  return std::unique_ptr<Method>{func(size, data)};
}


void method_serialization_deleter(MethodSerial *p) {
  hpx_free_registered(p);
}


/////////////////////////////////////////////////////////////////////
// In the public interface to dashmm
/////////////////////////////////////////////////////////////////////


bool register_method(int type, hpx_action_t creator) {
  int nlocs = hpx_get_num_ranks();
  hpx_addr_t checker = hpx_lco_reduce_new(nlocs, sizeof(int),
                                          int_sum_ident_op, int_sum_op);
  assert(checker != HPX_NULL);
  hpx_bcast_lsync(register_method_action, HPX_NULL, &type, &creator, &checker);
  int checkval{0};
  hpx_lco_get(checker, sizeof(int), &checkval);
  return (checkval == 0);
}


MethodSerialPtr method_serialization_allocator(size_t size) {
  MethodSerial *p = static_cast<MethodSerial *>(hpx_malloc_registered(size));
  return MethodSerialPtr{p, method_serialization_deleter};
}


} // namespace dashmm
