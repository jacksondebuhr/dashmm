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


struct ExpansionTableRow {
  hpx_action_t create;
  hpx_action_t interpret;
}


//map from type to description of that type
std::map<int, ExpansionTableRow> expansion_table_;


//parameter types for user-marshalled actions


/////////////////////////////////////////////////////////////////////
// Actions
/////////////////////////////////////////////////////////////////////


int register_expansion_handler(int type, hpx_action_t creator,
                               hpx_action_t interpreter, hpx_addr_t check) {
  int checkval = 0;
  if (expansion_table_.count(type) != 0) {
    int checkval = 1;
  } else {
    expansion_table_[type] = ExpansionTableRow{creator, interpreter};
  }
  hpx_lco_set_lsync(check, sizeof(int), &checkval, HPX_NULL);
}
HPX_ACTION(HPX_DEFAULT, 0,
           register_expansion_action, register_expansion_handler,
           HPX_INT, HPX_ACTION_T, HPX_ACTION_T, HPX_ADDR);


/////////////////////////////////////////////////////////////////////
// Internal Routines
/////////////////////////////////////////////////////////////////////


std::unique_ptr<Expansion> interpret_expansion(int type, void *data,
                                               size_t size) {
  auto entry = expansion_table_.find(type);
  if (entry == expansion_table_.end()) {
    return std::unique_ptr<Expansion>{nullptr};
  }
  expansion_creation_function_t func =
      reinterpret_cast<expansion_creation_function_t>(
        hpx_action_get_handler(entry->second.interpret)
      );
  return std::unique_ptr<Expansion>{func(size, data)};
}

std::unique_ptr<Expansion> create_expansion(int type, Point center) {
  auto entry = expansion_table_.find(type);
  if (entry == expansion_table_.end()) {
    return nullptr;
  }
  expansion_creation_function_t func =
      reinterpret_cast<expansion_creation_function_t>(
        hpx_action_get_handler(entry->second.create)
      );
  return std::unique_ptr<Expansion>{func(center)};
}


/////////////////////////////////////////////////////////////////////
// Public Interface
/////////////////////////////////////////////////////////////////////


//TODO: extend this to have both in-HPX and outside-HPX versions.
int register_expansion(int type, hpx_action_t creator,
                                 hpx_action_t interpreter) {
  if (type < kFirstUserExpansionType || type > kLastUserExpansionType) {
    return false;
  }

  int nlocs = hpx_get_num_ranks();
  hpx_addr_t checker = hpx_lco_reduce_new(nlocs, sizeof(int),
                                          int_sum_ident_op, int_sum_op);
  assert(checker != HPX_NULL);
  hpx_bcast_lsync(register_expansion_action, HPX_NULL, &type, &creator,
                  &interpreter, &checker);
  int checkval{0};
  hpx_lco_get(checker, sizeof(int), &checkval);
  return (checkval == 0);
}


} // namespace dashmm
