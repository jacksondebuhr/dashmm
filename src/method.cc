#include "method.h"

#include <map>

#include "hpx/hpx.h"
#include "libsync/sync.h"

//other dashmm stuff


namespace dashmm {


constexpr int kFirstMethodType = 0;
constexpr int kLastMethodType = 999;


//This is a global quantity containing the next available type for methods
// this is an integer. This is only used from locality zero. Other localities
// will send requests to locality zero for a new type identifier.
int next_available_method_ = kFirstMethodType;


//The mapping from type to table entries
std::map<int, MethodDesc> method_table_;


//Types for user-marshalled actions
struct register_method_params_t {
  int type;
  MethodDesc desc;
};


/////////////////////////////////////////////////////////////////////
// Actions
/////////////////////////////////////////////////////////////////////


int register_method_handler(register_method_params_t *parms, size_t size) {
  assert(method_table_[parms->type].count() == 0
            && "Registering method over existing");
  method_table_[parms->type] = parms->desc;
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
           register_method_action, register_method_handler,
           HPX_POINTER, HPX_SIZE_T);


int request_next_method_identifier_handler(void *UNUSED) {
  int retval = sync_fadd(&next_available_method_, 1, SYNC_ACQ_REL);
  assert(retval <= kLastMethodType && "Out of methods, somehow...");
  HPX_THREAD_CONTINUE(retval);
}
HPX_ACTION(HPX_DEFAULT, 0, request_next_method_identifier_action
           request_next_method_identifier_handler, HPX_POINTER);


int method_generate_handler(int expand_type, hpx_action_t genfunc) {
  hpx_addr_t source_node = hpx_thread_current_target();
  SourceNode node{source_node};
  Expansion expand{expand_type, Point{0.0, 0.0, 0.0}, false};
  generate_handler_t func = hpx_action_get_handler(genfunc);
  func(node, expand);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           method_generate_action, method_generate_handler,
           HPX_INT, HPX_ACTION_T);


int method_aggregate_handler(int expand_type, hpx_action_t aggfunc) {
  hpx_addr_t source_node = hpx_thread_current_target();
  SourceNode node{source_node};
  Expansion expand{expand_type, Point{0.0, 0.0, 0.0}, false};
  aggregate_handler_t func = hpx_action_get_handler(aggfunc);
  func(node, expand);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           method_aggregate_action, method_aggregate_handler,
           HPX_INT, HPX_ACTION_T);


int method_inherit_handler(int expand_type, hpx_action_t inhfunc,
                           int which_child) {
  hpx_addr_t target_node = hpx_thread_current_target();
  TargetNode node{target_node};
  Expansion expand{expand_type, Point{0.0, 0.0, 0.0}, false};
  inherit_handler_t func = hpx_action_get_handler(inhfunc);
  func(node, expand, which_child);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           method_inherit_action, method_inherit_handler,
           HPX_INT, HPX_ACTION_T, HPX_INT);


/////////////////////////////////////////////////////////////////////
// Internal routines
/////////////////////////////////////////////////////////////////////


const MethodDesc &get_method_desc(int type) {
  auto entry = method_table_.find(type);
  if (entry == method_table_.end()) {
    return method_table_[kFirstMethodType];
  } else {
    return *entry;
  }
}


int request_next_method_identifier() {
  int retval{0};
  hpx_call_sync(HPX_THERE(0), request_next_method_identifier_action,
                &retval, sizeof(retval), nullptr);
  return retval;
}


/////////////////////////////////////////////////////////////////////
// In the public interface to dashmm
/////////////////////////////////////////////////////////////////////


Method::Method(int type, std::vector<double> parms) {
  type_ = type;
  table_ = get_method_desc(type);
  if (table_.compatible_with_function == HPX_ACTION_NULL) {
    type_ = kFirstMethodType;
  } else {
    params_ = parms;
  }
}


bool Method::compatible_with(const Expansion &expand) const {
  compatible_with_handler_t func =
            hpx_action_get_handler(table_.compatible_with_function);
  return func(expand);
}


hpx_addr_t Method::generate(hpx_addr_t sync, SourceNode &curr,
                            const Expansion &expand) const {
  hpx_addr_t retval{sync};
  if (sync == HPX_NULL) {
    retval = hpx_lco_future_new(0);
    assert(retval != HPX_NULL);
  }

  int type{expand.type()};
  hpx_action_t genfunc{table_.generate_function};
  hpx_call(curr.data(), method_generate_action, retval, &type, &genfunc);

  return retval;
}


hpx_addr_t Method::aggregate(hpx_addr_t sync, SourceNode &curr,
                     const Expansion &expand) const {
  hpx_addr_t retval{sync};
  if (sync == HPX_NULL) {
    retval = hpx_lco_future_new(0);
    assert(retval != HPX_NULL);
  }

  int type{expand.type()};
  hpx_action_t aggfunc{table_.aggregate_function};
  hpx_call(curr.data(), method_aggregate_action, retval, &type, &aggfunc);

  return retval;
}


hpx_addr_t Method::inherit(hpx_addr_t sync, TargetNode &curr,
                           const Expansion &expand, size_t which_child) const {
  hpx_action_t retval{sync};
  if (sync == HPX_NULL) {
    retval = hpx_lco_future_new(0);
    assert(retval != HPX_NULL);
  }

  int type{expand.type()};
  hpx_action_t inhfunc{table_.inherit_function};
  hpx_call(curr.data(), method_inherit_action, retval, &type, &inhfunc,
           &which_child);

  return retval;
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


int register_method(MethodDesc desc) {
  register_method_params_t parms{};
  parms.type = request_next_method_identifier();
  parms.desc = desc;
  hpx_bcast_rsync(register_method_action, &parms, sizeof(parms));
  return parms.type;
}


} // namespace dashmm
