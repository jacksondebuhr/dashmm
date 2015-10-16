#include "expansion.h"

#include <map>

#include <hpx/hpx.h>

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


int register_expansion_handler(register_expansion_params_t *parms,
                               size_t size) {
  assert(expansion_table_[parms->type].count() == 0
            && "Registering expansion over existing");
  parms->desc.core_data = nullptr;
  if (parms->desc.core_size) {
    coregen_function_t func = hpx_action_gen_handler(parms->coregen);
    parms->desc.core_data = func(parms->desc.size);
  }
  expansion_table_[parms->type] = parms->desc;
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
           register_expansion_action, register_expansion_handler,
           HPX_POINTER, HPX_SIZE_T);


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


ExpansionRef globalize_expansion(Expansion *exp, hpx_addr_t where) {
  if (exp == nullptr) {
    return ExpansionRef{HPX_NULL};
  }
  ExpansionSerialPtr serial = met->serialize();
  size_t size = serial->size + sizeof(ExpansionSerial);
  hpx_addr_t data = hpx_gas_alloc_local_at_sync(size, 0, where);
  assert(data != HPX_NULL);
  hpx_gas_memput_rsync(data, serial.get(), size);
  return data;
}


} // namespace dashmm
