#include "expansion.h"

#include <map>

#include <hpx/hpx.h>

#include "include/expansionref.h"
#include "include/reductionops.h"


/// \file src/expansion.cc
/// \brief Implementation of expansion general framework


namespace dashmm {


/// The lowest allowed expansion identifier for internal use
constexpr int kFirstExpansionType = 2000;

/// The highest allowed expansion identifier for internal use
constexpr int kLastExpansionType = 2999;

/// The first allowed expansion identifier for user defined expansions
constexpr int kFirstUserExpansionType = 3000;

/// the highest allowed expansion identifier for user defined expansions
constexpr int kLastUserExpansionType = 3999;


constexpr int kAllocateExpansionTable = 0;
constexpr int kDeleteExpansionTable = 1;


/// A row of the table serving expansion creation functions
struct ExpansionTableRow {
  hpx_action_t create;
  hpx_action_t interpret;
}


/// A mapping from expansion type to creation functions
///
/// This is a local object relative to HPX-5. What this means is that each
/// (traditional) process will have a copy of this table. Registering an
/// action amounts to a broadcast.
std::map<int, ExpansionTableRow> *expansion_table_;


/////////////////////////////////////////////////////////////////////
// Actions
/////////////////////////////////////////////////////////////////////


/// Action to register an expansion with DASHMM
///
/// This action is invoked with a broadcast.
///
/// \param type - the type identifier being registered
/// \param creator - the HPX function implementing creation of that type
/// \param interpreter - the HPX function implementing interpretation of that
///                      type
/// \param check - the address of a reduction LCO that will be used for error
///                reporting. It will have a nonzero contribution on error.
///
/// \returns - HPX_SUCCESS
int register_expansion_handler(int type, hpx_action_t creator,
                               hpx_action_t interpreter, hpx_addr_t check) {
  int checkval = 0;
  if (expansion_table_->count(type) != 0) {
    int checkval = 1;
  } else {
    (*expansion_table_)[type] = ExpansionTableRow{creator, interpreter};
  }
  hpx_lco_set_lsync(check, sizeof(int), &checkval, HPX_NULL);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           register_expansion_action, register_expansion_handler,
           HPX_INT, HPX_ACTION_T, HPX_ACTION_T, HPX_ADDR);


int manage_expansion_table_handler(int operation) {
  if (operation == kAllocateExpansionTable) {
    expansion_table_ = new std::map<int, ExpansionTableRow>{};
  } else if (operation == kDeleteExpansionTable) {
    delete expansion_table_;
  }
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           manage_expansion_table_action, manage_expansion_table_handler,
           HPX_INT);


/////////////////////////////////////////////////////////////////////
// Internal Routines
/////////////////////////////////////////////////////////////////////


std::unique_ptr<Expansion> interpret_expansion(int type, void *data,
                                               size_t size) {
  auto entry = expansion_table_->find(type);
  if (entry == expansion_table_->end()) {
    return std::unique_ptr<Expansion>{nullptr};
  }
  expansion_creation_function_t func =
      reinterpret_cast<expansion_creation_function_t>(
        hpx_action_get_handler(entry->second.interpret)
      );
  return std::unique_ptr<Expansion>{func(size, data)};
}


std::unique_ptr<Expansion> create_expansion(int type, Point center) {
  auto entry = expansion_table_->find(type);
  if (entry == expansion_table_->end()) {
    return nullptr;
  }
  expansion_creation_function_t func =
      reinterpret_cast<expansion_creation_function_t>(
        hpx_action_get_handler(entry->second.create)
      );
  return std::unique_ptr<Expansion>{func(center)};
}


void init_expansion_table() {
  int input = kAllocateExpansionTable;
  hpx_bcast_rsync(manage_expansion_table_handler, &input);
}


void fini_expansion_table() {
  int input = kAllocateExpansionTable;
  hpx_bcast_rsync(manage_expansion_table_handler, &input);
}


/////////////////////////////////////////////////////////////////////
// Public Interface
/////////////////////////////////////////////////////////////////////


ReturnCode register_expansion(int type, hpx_action_t creator,
                                 hpx_action_t interpreter) {
  if (type < kFirstUserExpansionType || type > kLastUserExpansionType) {
    return kDomainError;
  }

  int nlocs = hpx_get_num_ranks();
  hpx_addr_t checker = hpx_lco_reduce_new(nlocs, sizeof(int),
                                          int_sum_ident_op, int_sum_op);
  assert(checker != HPX_NULL);
  hpx_bcast_lsync(register_expansion_action, HPX_NULL, &type, &creator,
                  &interpreter, &checker);
  int checkval{0};
  hpx_lco_get(checker, sizeof(int), &checkval);
  return (checkval == 0 ? kSuccess : kDomainError);
}


} // namespace dashmm
