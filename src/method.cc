#include "method.h"

#include <map>

#include "hpx/hpx.h"
#include "libsync/sync.h"

//other dashmm stuff


namespace dashmm {


constexpr int first_method_type_ = 0;
constexpr int last_method_type_ = 999;


//This is a global quantity containing the next available type for methods
// this is an integer. This is only used from locality zero. Other localities
// will send requests to locality zero for a new type identifier.
int next_available_method_ = first_method_type_;


//The mapping from type to table entries
std::map<int, MethodDesc> method_table_;


/////////////////////////////////////////////////////////////////////
// Actions
/////////////////////////////////////////////////////////////////////


int register_method_handler(int type, size_t params, hpx_action_t compat,
                            hpx_action_t gen, hpx_action_t agg,
                            hpx_action_t inherit, hpx_action_t proc,
                            hpx_action_t reftest) {
  MethodDesc desc{params, compat, gen, agg, inherit, proc, reftest};
  assert(method_table[type].count() == 0 && "Registering method over existing");
  method_table_[type] = desc;
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, register_method_action, register_method_handler,
           HPX_INT, HPX_SIZE_T, HPX_ACTION_T, HPX_ACTION_T, HPX_ACTION_T,
           HPX_ACTION_T, HPX_ACTION_T, HPX_ACTION_T);


int request_next_method_identifier_handler(void *UNUSED, size_t UNWANTED) {
  int retval = sync_fadd(&next_available_method_, 1, SYNC_ACQ_REL);
  assert(retval <= last_method_type_ && "Out of methods, somehow...");
  HPX_THREAD_CONTINUE(retval);
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED, request_next_method_identifier_action
           request_next_method_identifier_handler, HPX_POINTER, HPX_SIZE_T);


/////////////////////////////////////////////////////////////////////
// Internal routines
/////////////////////////////////////////////////////////////////////


const MethodDesc &get_method_desc(int type) {
  auto entry = method_table_.find(type);
  if (entry == method_table_.end()) {
    return method_table_[first_method_type_];
  } else {
    return *entry;
  }
}


int request_next_method_identifier() {
  int retval{0};
  hpx_call_sync(HPX_THERE(0), request_next_method_identifier_action,
                &retval, sizeof(retval), nullptr, 0);
  return retval;
}


/////////////////////////////////////////////////////////////////////
// In the public interface to dashmm
/////////////////////////////////////////////////////////////////////


int register_method(size_t params, hpx_action_t compat,
                    hpx_action_t gen, hpx_action_t agg, hpx_action_t inherit,
                    hpx_action_t proc, hpx_action_t reftest) {
  int retval = request_next_method_identifier();
  hpx_bcast_rsync(register_method_action, &retval, &params, &compat, &gen,
                  &agg, &inherit, &proc, &reftest);
  return retval;
}


Method::Method(int type, std::vector<double> parms) {
  type_ = type;
  table_ = get_method_desc(type);
  if (table_.compatible_with_function == HPX_ACTION_NULL) {
    type_ = first_method_type_;
  } else {
    params_ = parms;
  }
}


bool Method::valid() const {
  return type_ != first_method_type_;
}


hpx_addr_t Method::compatible_with(hpx_addr_t sync,
                                   const Expansion *expand) const {
  hpx_addr_t retval = sync;
  if (retval == HPX_NULL) {
    retval = hpx_lco_future_new(sizeof(bool));
  }

  hpx_call();

  return retval;
}


hpx_addr_t Method::generate(hpx_addr_t sync, SourceNode *curr,
                    const Expansion *expand) const {
  //
}


hpx_addr_t Method::aggregate(hpx_addr_t sync, SourceNode *curr,
                     const Expansion *expand) const {
  //
}


hpx_addr_t Method::inherit(hpx_addr_t sync, TargetNode *curr,
                           const Expansion *expand, size_t which_child) const {
  //
}


hpx_addr_t Method::process(hpx_addr_t sync, TargetNode *curr,
                           std::vector<SourceNode *> &consider,
                           bool curr_is_leaf) const {
  //
}


hpx_addr_t Method::refine_test(hpx_addr_t sync, bool same_sources_and_targets,
                            const TargetNode *curr,
                            const std::vector<SourceNode *> &consider) const {
  //
}


} // namespace dashmm
