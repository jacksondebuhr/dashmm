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


//The Method's data
struct MethodData {
  int type;
  const MethodDesc &table;
  size_t param_count;
  double params[];
};


//Types for user-marshalled actions
struct register_method_params_t {
  int type;
  MethodDesc desc;
};

struct method_process_params_t {
  hpx_action_t procfunc;
  bool curr_is_leaf;
  size_t consider_size;
  hpx_addr_t consider[];
};

struct method_refine_test_params_t {
  hpx_action_t refine_func;
  bool same_sources_and_targets;
  size_t consider_size;
  hpx_addr_t consider[];
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


int method_process_handler(method_process_params_t *parms, size_t size) {
  hpx_addr_t target_node = hpx_thread_current_target();
  TargetNode node{target_node};
  std::vector<SourceNode> consider(parms->consider_size);
  for (size_t i = 0; i < parms->consider_size; ++i) {
    consider[i] = SourceNode{parms->consider[i]};
  }
  process_handler_t func = hpx_action_get_handler(parms->procfunc);
  func(node, consider, parms->curr_is_leaf);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
           method_process_action, method_process_handler,
           HPX_POINTER, HPX_SIZE_T);


int method_refine_test_handler(method_refine_test_params_t *parms,
                               size_t size) {
  hpx_addr_t target_node = hpx_thread_current_target();
  TargetNode node{target_node};
  std::vector<SourceNode> consider(parms->consider_size);
  for (size_t i = 0; i < parms.consider_size; ++i) {
    consider[i] = SourceNode{parms->consider[i]};
  }
  refine_test_handler_t func = parms->refine_func;
  bool result = func(parms->same_sources_and_targets, node, consider);
  HPX_THREAD_CONTINUE(result);
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
           method_refine_test_action, method_refine_test_handler,
           HPX_POINTER, HPX_SIZE_T);


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


size_t size_of_method_process_params(size_t count) {
  return sizeof(method_process_params_t) + sizeof(hpx_addr_t) * count;
}


size_t size_of_method_refine_test_params(size_t count) {
  return sizeof(method_refine_test_params_t) + sizeof(hpx_addr_t) * count;
}


size_t size_of_method_data(MethodData *mdata) {
  return sizeof(MethodData) + mdata->param_count * sizeof(double);
}


MethodData *method_data_pin(hpx_addr_t data) {
  MethodData *local{nullptr};
  assert(hpx_gas_try_pin(data, (void **)&local));
  return local;
}


void method_data_unpin(hpx_addr_t data) {
  hpx_gas_unpin(data);
}


void method_is_valid(MethodData *method) {
  return method->type != kFirstMethodType;
}


/////////////////////////////////////////////////////////////////////
// In the public interface to dashmm
/////////////////////////////////////////////////////////////////////


Method::Method(int type, std::vector<double> parms) {
  data_ = hpx_gas_alloc(1, size_of_method_data(parms.size()), 0,
                        HPX_DIST_TYPE_LOCAL);
  assert(data_ != HPX_NULL);

  MethodData *method = method_data_pin(data_);
  method->type = type;
  //NOTE: This only works in SMP. In distributed, we will need to keep refinding
  // this reference when we use the MethodData
  method->table = get_method_desc(type);
  if (method->table.compatible_with_function == HPX_ACTION_NULL) {
    method->type = kFirstMethodType;
  } else {
    method->param_count = parms.size();
    for (size_t i = 0; i < parms.size(); ++i) {
      method->params[i] = parms[i];
    }
  }
  method_data_unpin(data_);
}


bool Method::valid() const {
  return type() != kFirstMethodType;
}


int Method::type() const {
  MethodData *method = method_data_pin(data_);
  int retval{method->type};
  method_data_unpin(data_);
  return retval;
}


bool Method::compatible_with(const Expansion &expand) const {
  MethodData *method = method_data_pin(data_);
  assert(method_is_valid(method));
  compatible_with_handler_t func =
            hpx_action_get_handler(method->table_.compatible_with_function);
  method_data_unpin(data_);
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
  MethodData *method = method_data_pin(data_);
  hpx_action_t genfunc{method->table_.generate_function};
  method_data_unpin(data_);
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
  MethodData *method = method_data_pin(data_);
  hpx_action_t aggfunc{method->table_.aggregate_function};
  method_data_unpin(data_);
  hpx_call(curr.data(), method_aggregate_action, retval, &type, &aggfunc);

  return retval;
}


hpx_addr_t Method::inherit(hpx_addr_t sync, TargetNode &curr,
                           const Expansion &expand, size_t which_child) const {
  hpx_addr_t retval{sync};
  if (sync == HPX_NULL) {
    retval = hpx_lco_future_new(0);
    assert(retval != HPX_NULL);
  }

  int type{expand.type()};
  MethodData *method = method_data_pin(data_);
  hpx_action_t inhfunc{method->table_.inherit_function};
  method_data_unpin(data_);
  hpx_call(curr.data(), method_inherit_action, retval, &type, &inhfunc,
           &which_child);

  return retval;
}


hpx_addr_t Method::process(hpx_addr_t sync, TargetNode &curr,
                           std::vector<SourceNode> &consider,
                           bool curr_is_leaf) const {
  hpx_addr_t retval{sync};
  if (sync == HPX_NULL) {
    retval = hpx_lco_future_new(0);
    assert(retval != HPX_NULL);
  }

  size_t parms_size = size_of_method_process_params(consider.size());
  method_process_params_t *parms = malloc(parms_size);
  assert(parms != nullptr);
  parms->proc_func = table_.process_function;
  parms->curr_is_leaf = curr_is_leaf;
  parms->consider_size = consider.size();
  for (size_t i = 0; i < consider.size(); ++i) {
    parms->consider[i] = consider[i].data();
  }
  hpx_call(curr.data(), method_process_action, retval, parms, parms_size);
  free(parms);

  return retval;
}


hpx_addr_t Method::refine_test(hpx_addr_t sync, bool same_sources_and_targets,
                            const TargetNode &curr,
                            const std::vector<SourceNode> &consider) const {
  hpx_addr_t retval{sync};
  if (sync == HPX_NULL) {
    retval = hpx_lco_future_new(sizeof(int));
    assert(retval != HPX_NULL);
  }

  size_t parms_size = size_of_method_refine_test_params(consider.size());
  method_refine_test_params_t *parms = malloc(parms_size);
  assert(parms != nullptr);
  parms->refine_func = table_.refine_test_function;
  parms->same_sources_and_targets = same_sources_and_targets;
  parms->consider_size = consider.size();
  for (size_t i = 0; i < consider.size(); ++i) {
    parms->consider[i] = consider[i].data();
  }
  hpx_call(curr.data(), method_refine_test_action, retval, parms, parms_size);
  free(parms);

  return retval;
}


int register_method(MethodDesc desc) {
  register_method_params_t parms{};
  parms.type = request_next_method_identifier();
  parms.desc = desc;
  hpx_bcast_rsync(register_method_action, &parms, sizeof(parms));
  return parms.type;
}


} // namespace dashmm
