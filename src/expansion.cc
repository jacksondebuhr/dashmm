#include "expansion.h"

#include <map>

#include <hpx/hpx.h>

#include "include/expansionref.h"
#include "include/reductionops.h"


namespace dashmm {


//The range of Expansion type identifiers
constexpr int kFirstExpansionType = 1000;
constexpr int kLastExpansionType = 1999;
constexpr int kFirstUserExpansionType = 2000;
constexpr int kLastUserExpansionType = 2999;


//map from type to description of that type
std::map<int, hpx_action_t> expansion_table_;


//parameter types for user-marshalled actions


/////////////////////////////////////////////////////////////////////
// Actions
/////////////////////////////////////////////////////////////////////


int register_expansion_handler(int type, hpx_action_t creator,
                               hpx_addr_t check) {
  int checkval = 0;
  if (expansion_table_.count(type) != 0) {
    int checkval = 1;
  } else {
    expansion_table_[type] = creator;
  }
  hpx_lco_set_lsync(check, sizeof(int), &checkval, HPX_NULL);
}
HPX_ACTION(HPX_DEFAULT, 0,
           register_expansion_action, register_expansion_handler,
           HPX_INT, HPX_ACTION_T, HPX_ADDR);


/////////////////////////////////////////////////////////////////////
// Internal Routines
/////////////////////////////////////////////////////////////////////


std::unique_ptr<Expansion> create_expansion(int type, size_t size, void *data) {
  auto entry = expansion_table_.find(type);
  if (entry == expansion_table_.end()) {
    return std::unique_ptr<Expansion>{nullptr};
  }
  expansion_creation_function_t func =
      reinterpret_cast<expansion_creation_function_t>(
        hpx_action_get_handler(entry->second)
      );
  return std::unique_ptr<Expansion>{func(size, data)};
}


void expansion_serialization_deleter(ExpansionSerial *p) {
  hpx_free_registered(p);
}


/////////////////////////////////////////////////////////////////////
// Public Interface
/////////////////////////////////////////////////////////////////////


int register_expansion(int type, hpx_action_t creator) {
  int nlocs = hpx_get_num_ranks();
  hpx_addr_t checker = hpx_lco_reduce_new(nlocs, sizeof(int),
                                          int_sum_ident_op, int_sum_op);
  assert(checker != HPX_NULL);
  hpx_bcast_lsync(register_expansion_action, HPX_NULL, &type, &creator,
                  &checker);
  int checkval{0};
  hpx_lco_get(checker, sizeof(int), &checkval);
  return (checkval == 0);
}


ExpansionSerialPtr expansion_serialization_allocator(size_t size) {
  ExpansionSerial *p = static_cast<ExpansionSerial *>(
      hpx_malloc_registered(size));
  return ExpansionSerialPtr{p, expansion_serialization_deleter};
}


} // namespace dashmm
