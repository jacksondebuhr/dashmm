/// \file src/node.cc
/// \brief implementation of SourceNode and TargetNode


#include "include/node.h"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>

#include <algorithm>

#include <hpx/hpx.h>

#include "include/array.h"
#include "include/domaingeometry.h"
#include "include/index.h"
#include "include/method.h"
#include "include/methodref.h"


namespace dashmm {


/// The data for source nodes.
///
/// The data stored for SourceNode objects. This will be saved in a block in
/// the GAS.
struct SourceNodeData {
  /// The geometry of the root for the tree of which this node is a part.
  DomainGeometry root_geo;
  /// The index giving which subdivision of the root this node is.
  Index idx;
  /// The global address of the parent of this node.
  hpx_addr_t parent;
  /// The global addresses of the children of this node.
  hpx_addr_t child[8];
  /// The global address of the expansion object for this node.
  ExpansionRef expansion;
  /// The global address of the method object for this node.
  hpx_addr_t method;
  /// A reference to the sources for this node. If this is an internal node,
  /// this will be an invalid SourceRef.
  SourceRef sources;
};


/// The data for target nodes.
///
/// The data stored for TargetNode objects. This will be saved in a block in
/// the GAS.
struct TargetNodeData {
  /// The geometry of the root for the tree of which this node is a part.
  DomainGeometry root_geo;
  /// The index giving which subdivision of the root this node is.
  Index idx;
  /// The global address of the parent of this node.
  hpx_addr_t parent;
  /// The global addresses of the children of this node.
  hpx_addr_t child[8];
  /// The global address of the expansion object for this node.
  ExpansionRef expansion;
  /// The global address of the method object for this node.
  MethodRef method;
  /// A reference to the targets for this node.
  TargetRef targets;
};


/////////////////////////////////////////////////////////////////////
// Actions and other HPX-5 related
/////////////////////////////////////////////////////////////////////


int node_delete_self_handler(hpx_addr_t gate) {
  if (gate != HPX_NULL) {
    hpx_lco_delete_sync(gate);
  }
  hpx_addr_t target = hpx_thread_current_target();
  hpx_gas_free_sync(target);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0, node_delete_self_action, node_delete_self_handler,
           HPX_ADDR);


int source_node_delete_handler(hpx_addr_t data);
HPX_ACTION(HPX_DEFAULT, 0,
           source_node_delete_action, source_node_delete_handler,
           HPX_ADDR);

int source_node_delete_handler(hpx_addr_t data) {
  SourceNodeData *local{nullptr};
  assert(hpx_gas_try_pin(data, (void **)&local));

  local->expansion.destroy();
  local->sources.destroy();

  int count{0};
  for (int i = 0; i < 8; ++i) {
    if (local->child[i] != HPX_NULL) {
      ++count;
    }
  }

  hpx_addr_t done{HPX_NULL};
  if (count) {
    done = hpx_lco_and_new(count);
    for (int i = 0; i < 8; ++i) {
      if (local->child[i] != HPX_NULL) {
        hpx_call(local->child[i], source_node_delete_action, done,
                 &local->child[i]);
      }
    }
  }

  hpx_gas_unpin(data);

  //call when with current continuation
  // This relies on the fact that HPX_NULL as the gate means this will
  // be equivalent to hpx_call_cc.
  hpx_call_when_cc(done, data, node_delete_self_action, nullptr, nullptr,
                   &done);
}


int target_node_delete_handler(hpx_addr_t data);
HPX_ACTION(HPX_DEFAULT, 0,
           target_node_delete_action, target_node_delete_handler,
           HPX_ADDR);

int target_node_delete_handler(hpx_addr_t data) {
  TargetNodeData *local{nullptr};
  assert(hpx_gas_try_pin(data, (void **)&local));

  local->expansion.destroy();
  local->targets.destroy();

  int count{0};
  for (int i = 0; i < 8; ++i) {
    if (local->child[i] != HPX_NULL) {
      ++count;
    }
  }

  hpx_addr_t done{HPX_NULL};
  if (count) {
    done = hpx_lco_and_new(count);
    for (int i = 0; i < 8; ++i) {
      if (local->child[i] != HPX_NULL) {
        hpx_call(local->child[i], target_node_delete_action, done,
                 &local->child[i]);
      }
    }
  }

  hpx_gas_unpin(data);

  //call when with current continuation
  // This relies on the fact that HPX_NULL as the gate means this will
  // be equivalent to hpx_call_cc.
  hpx_call_when_cc(done, data, node_delete_self_action, nullptr, nullptr,
                   &done);
}


int source_node_child_partition_done_handler(SourceNodeData *node,
                                              hpx_addr_t partdone,
                                              int type,
                                              hpx_addr_t expand) {
  ExpansionRef expref{type, expand};
  MethodRef method{node->method};
  SourceNode snode{hpx_thread_current_target()};

  method.aggregate(snode, expref);

  node->expansion.finalize();

  hpx_lco_delete_sync(partdone);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           source_node_child_partition_done_action,
           source_node_child_partition_done_handler,
           HPX_POINTER, HPX_ADDR, HPX_INT, HPX_ADDR);


struct SourceNodePartitionParams {
  hpx_addr_t partdone;
  int limit;
  int type;
  hpx_addr_t expand;
  int n_sources;
  Source sources[];
};

int source_node_partition_handler(SourceNodeData *node,
                                  SourceNodePartitionParams *parms,
                                  size_t bytes);
HPX_ACTION(HPX_DEFAULT, HPX_PINNED | HPX_MARSHALLED,
           source_node_partition_action, source_node_partition_handler,
           HPX_POINTER, HPX_POINTER, HPX_SIZE_T);

int source_node_partition_handler(SourceNodeData *node,
                                  SourceNodePartitionParams *parms,
                                  size_t bytes) {
  if (parms->n_sources <= parms->limit) {
    assert(node->sources.data() == HPX_NULL);
    node->sources = SourceRef(parms->sources, parms->n_sources);

    MethodRef method{node->method};
    SourceNode curr{hpx_thread_current_target()};
    ExpansionRef expref{parms->type, parms->expand}; //this is the prototype

    //This will cause the creation of a new expansion, that will be globalized
    // and then set as the expansion for this node. It will be created with
    // initial data, and so there is no need to schedule any contributions to
    // the expansion.
    method.generate(curr, expref);

    //We are done scheduling contributions to the expansion for this node.
    node->expansion.finalize();

    hpx_lco_set(parms->partdone, 0, nullptr, HPX_NULL, HPX_NULL);
    return HPX_SUCCESS;
  }

  //NOTE: We make use of the SMP only form here. This will not work in
  // distributed.

  //partition sources
  Source *splits[9]{};
  splits[0] = parms->sources;
  splits[8] = &parms->sources[parms->n_sources];

  Point cen{node->root_geo.center_from_index(node->idx)};
  double z_center = cen.z();
  auto z_comp = [&z_center](Source &a) {
    return a.position.z() < z_center;
  };
  splits[4] = std::partition(splits[0], splits[8], z_comp);

  double y_center = cen.y();
  auto y_comp = [&y_center](Source &a) {
    return a.position.y() < y_center;
  };
  splits[2] = std::partition(splits[0], splits[4], y_comp);
  splits[6] = std::partition(splits[4], splits[8], y_comp);

  double x_center = cen.x();
  auto x_comp = [&x_center](Source &a) {
    return a.position.x() < x_center;
  };
  splits[1] = std::partition(splits[0], splits[2], x_comp);
  splits[3] = std::partition(splits[2], splits[4], x_comp);
  splits[5] = std::partition(splits[4], splits[6], x_comp);
  splits[7] = std::partition(splits[6], splits[8], x_comp);

  //Find some counts
  Source *cparts[8]{};
  int n_per_child[8]{0, 0, 0, 0, 0, 0, 0, 0};
  int n_children{0};
  int n_offset{0};
  for (int i = 0; i < 8; ++i) {
    n_per_child[i] = splits[i + 1] - splits[i];
    if (n_per_child[i]) {
      ++n_children;
      cparts[i] = &parms->sources[n_offset];
    } else {
      cparts[i] = nullptr;
    }
    n_offset += n_per_child[i];
  }

  hpx_addr_t childpartdone = hpx_lco_and_new(n_children);
  assert(childpartdone != HPX_NULL);

  for (int i = 0; i < 8; ++i) {
    if (n_per_child[i] == 0) {
      node->child[i] = HPX_NULL;
      continue;
    }

    Index cidx{node->idx.child(i)};
    SourceNode thisnode{hpx_thread_current_target()};
    SourceNode kid{node->root_geo, cidx,
                   node->method, &thisnode};
    node->child[i] = kid.data();

    size_t argsz = sizeof(SourceNodePartitionParams)
                   + sizeof(Source) * n_per_child[i];
    SourceNodePartitionParams *args =
        static_cast<SourceNodePartitionParams *>(malloc(argsz));
    assert(args != nullptr);
    args->partdone = childpartdone;
    args->limit = parms->limit;
    args->type = parms->type;
    args->expand = parms->expand;
    args->n_sources = n_per_child[i];
    memcpy(args->sources, cparts[i], sizeof(Source) * n_per_child[i]);
    hpx_call(kid.data(), source_node_partition_action, HPX_NULL, args, argsz);
    free(args);
  }

  hpx_call_when(childpartdone, hpx_thread_current_target(),
                source_node_child_partition_done_action, parms->partdone,
                &childpartdone, &parms->type, &parms->expand);

  return HPX_SUCCESS;
}


struct TargetNodePartitionParams {
  hpx_addr_t done;
  bool same_sources_and_targets;
  hpx_addr_t parts;
  int n_parts;
  int n_parts_total;
  int limit;
  ExpansionRef expansion;
  int which_child;
  int n_consider;
  hpx_addr_t consider[];
};

size_t target_node_partition_params_size(int n_consider) {
  return sizeof(TargetNodePartitionParams) + n_consider * sizeof(hpx_addr_t);
}

TargetNodePartitionParams *target_node_partition_params_alloc(int n_consider) {
  TargetNodePartitionParams *retval = static_cast<TargetNodePartitionParams *>(
      malloc(target_node_partition_params_size(n_consider)));
  if (retval) {
    retval->n_consider = n_consider;
  }
  return retval;
}

int target_node_partition_handler(TargetNodeData *node,
                                  TargetNodePartitionParams *parms,
                                  size_t bytes);
HPX_ACTION(HPX_DEFAULT, HPX_PINNED | HPX_MARSHALLED,
           target_node_partition_action, target_node_partition_handler,
           HPX_POINTER, HPX_POINTER, HPX_SIZE_T);

int target_node_partition_handler(TargetNodeData *node,
                                  TargetNodePartitionParams *parms,
                                  size_t bytes) {
  MethodRef method{node->method.data()};
  TargetNode curr{hpx_thread_current_target()};
  std::vector<SourceNode> consider{};
  for (int i = 0; i < parms->n_consider; ++i) {
    consider.push_back(SourceNode{parms->consider[i]});
  }

  bool refine = false;
  if (parms->n_parts > parms->limit) {
    refine = method.refine_test(parms->same_sources_and_targets, curr,
                                consider);
  }

  if (refine) {
    //partition
    Target *T{nullptr};
    assert(hpx_gas_try_pin(parms->parts, (void **)&T));
    Target *splits[9]{};
    splits[0] = T;
    splits[8] = &T[parms->n_parts];

    Point cen{node->root_geo.center_from_index(node->idx)};
    double z_center = cen.z();
    auto z_comp = [&z_center](Target &a) {
      return a.position.z() < z_center;
    };
    splits[4] = std::partition(splits[0], splits[8], z_comp);

    double y_center = cen.y();
    auto y_comp = [&y_center](Target &a) {
      return a.position.y() < y_center;
    };
    splits[2] = std::partition(splits[0], splits[4], y_comp);
    splits[6] = std::partition(splits[4], splits[8], y_comp);

    double x_center = cen.x();
    auto x_comp = [&x_center](Target &a) {
      return a.position.x() < x_center;
    };
    splits[1] = std::partition(splits[0], splits[2], x_comp);
    splits[3] = std::partition(splits[2], splits[4], x_comp);
    splits[5] = std::partition(splits[4], splits[6], x_comp);
    splits[7] = std::partition(splits[6], splits[8], x_comp);
    hpx_gas_unpin(parms->parts);

    hpx_addr_t cparts[8];
    int n_per_child[8]{0, 0, 0, 0, 0, 0, 0, 0};
    int n_children{0};
    int n_offset{0};
    for (int i = 0; i < 8; ++i) {
      n_per_child[i] = splits[i + 1] - splits[i];
      if (n_per_child[i]) {
        ++n_children;
        cparts[i] = hpx_addr_add(parms->parts, sizeof(Target) * n_offset,
                                        sizeof(Target) * parms->n_parts_total);
      } else {
        cparts[i] = HPX_NULL;
      }
      n_offset += n_per_child[i];
    }

    hpx_addr_t cdone = hpx_lco_and_new(n_children);
    assert(cdone != HPX_NULL);

    //set up the arguments to the partition actions; the constant parts
    TargetNodePartitionParams *args =
        target_node_partition_params_alloc(consider.size());
    size_t argssize = target_node_partition_params_size(consider.size());
    args->done = cdone;
    args->same_sources_and_targets = parms->same_sources_and_targets;
    args->n_parts_total = parms->n_parts_total;
    args->limit = parms->limit;
    args->expansion = parms->expansion;
    args->n_consider = consider.size();
    for (size_t i = 0; i < consider.size(); ++i) {
      args->consider[i] = consider[i].data();;
    }

    for (int i = 0; i < 8; ++i) {
      if (n_per_child[i] == 0) {
        node->child[i] = HPX_NULL;
        continue;
      }

      Index cidx{node->idx.child(i)};
      TargetNode thisnode{hpx_thread_current_target()};
      TargetNode kid{node->root_geo, cidx, node->method.data(), &thisnode};
      node->child[i] = kid.data();

      args->parts = cparts[i];
      args->n_parts = n_per_child[i];
      args->which_child = i;

      hpx_call(kid.data(), target_node_partition_action, HPX_NULL, args,
               argssize);
    }
    free(args);

    hpx_call_when(cdone, cdone, hpx_lco_delete_action, parms->done, nullptr, 0);
  } else {  // no refinement needed; so just set the input done LCO
    Target *targs{nullptr};
    assert(hpx_gas_try_pin(parms->parts, (void **)&targs));
    node->targets = TargetRef(targs, parms->n_parts);

    hpx_call_when(node->targets.data(), parms->done, hpx_lco_set_action,
                  HPX_NULL, nullptr, 0);
  }

  ExpansionRef expand{parms->expansion};
  method.inherit(curr, expand, parms->which_child);
  method.process(curr, consider, !refine);

  //At this point, all work on the current expansion will have been scheduled,
  // so we mark the LCO as such.
  expand.finalize();
  //Also, all the work on the targets for this node will have been scheduled as
  // well.
  if (!refine) {
    node->targets.finalize();
  }

  return HPX_SUCCESS;
}


int target_node_collect_results_handler(TargetNodeData *node,
                      hpx_addr_t user_array, size_t phi_offset);
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           target_node_collect_results_action,
           target_node_collect_results_handler,
           HPX_POINTER, HPX_ADDR, HPX_SIZE_T);

int target_node_collect_results_handler(TargetNodeData *node,
                      hpx_addr_t user_array, size_t phi_offset) {
  //count children
  int n_children{0};
  for (int i = 0; i < 8; ++i) {
    if (node->child[i] != HPX_NULL) {
      ++n_children;
    }
  }

  if (n_children == 0) {
    //pin targets of this node
    //HACK: This uses the size hack again.
    char *lcodata{nullptr};
    hpx_lco_getref(node->targets.data(), 1, (void **)&lcodata);
    Target *targets = reinterpret_cast<Target *>(lcodata + 4 * sizeof(int));

    //pin the user array
    ArrayMetaData *meta{nullptr};
    assert(hpx_gas_try_pin(user_array, (void **)&meta));
    char *user_data{nullptr};
    assert(hpx_gas_try_pin(meta->data, (void **)&user_data));

    //copy over
    for (int i = 0; i < node->targets.n(); ++i) {
      size_t idx = targets[i].index;
      char *record_base = &user_data[idx * meta->size];
      double *phi = reinterpret_cast<double *>(record_base + phi_offset);
      phi[0] = targets[i].phi.real();
      phi[1] = targets[i].phi.imag();
    }

    hpx_gas_unpin(meta->data);
    hpx_gas_unpin(user_array);
    hpx_lco_release(node->targets.data(), lcodata);

    return HPX_SUCCESS;
  } else {
    //internal, so set up calls for children
    hpx_addr_t coll_done = hpx_lco_and_new(n_children);
    assert(coll_done != HPX_NULL);

    for (int i = 0; i < 8; ++i) {
      if (node->child[i] != HPX_NULL) {
        hpx_call(node->child[i], target_node_collect_results_action, coll_done,
                 &user_array, &phi_offset);
      }
    }

    hpx_call_when_cc(coll_done, coll_done, hpx_lco_delete_action,
                     nullptr, nullptr, nullptr, 0);
  }
}


/////////////////////////////////////////////////////////////////////
// SourceNode implementation
/////////////////////////////////////////////////////////////////////


SourceNode::SourceNode(DomainGeometry g, Index idx,
                       hpx_addr_t method, SourceNode *parent) {
  data_ = hpx_gas_alloc_local(1, sizeof(SourceNodeData), 0);
  if (data_ == HPX_NULL) {
    return;
  }
  local_ = nullptr;

  pin();
  local_->root_geo = g;
  local_->idx = idx;
  if (parent) {
    local_->parent = parent->data();
  } else {
    local_->parent = HPX_NULL;
  }
  for (int i = 0; i < 8; ++i) {
    local_->child[i] = HPX_NULL;
  }
  local_->expansion = ExpansionRef{0, HPX_NULL};
  local_->method = method;
  local_->sources = SourceRef{};
}


SourceNode::~SourceNode() {
  unpin();
}


void SourceNode::destroy() {
  unpin();
  hpx_call_sync(data_, source_node_delete_action, nullptr, 0, &data_);
}


hpx_addr_t SourceNode::partition(Source *sources, int n_sources, int limit,
                                 int type, hpx_addr_t expand) {
  hpx_addr_t retval = hpx_lco_future_new(0);
  assert(retval != HPX_NULL);

  size_t bytes = sizeof(SourceNodePartitionParams) + n_sources * sizeof(Source);
  SourceNodePartitionParams *args =
      static_cast<SourceNodePartitionParams *>(malloc(bytes));
  assert(args);
  args->partdone = retval;
  args->limit = limit;
  args->type = type;
  args->expand = expand;
  args->n_sources = n_sources;
  memcpy(args->sources, sources, sizeof(Source) * n_sources);
  hpx_call(data_, source_node_partition_action, HPX_NULL,
           args, bytes);
  free(args);

  return retval;
}


bool SourceNode::is_leaf() const {
  pin();

  bool retval{true};
  for (int i = 0; i < 8; ++i) {
    if (local_->child[i] != HPX_NULL) {
      retval = false;
      break;
    }
  }

  return retval;
}


DomainGeometry SourceNode::root_geo() const {
  pin();
  return local_->root_geo;
}


Index SourceNode::index() const {
  pin();
  return local_->idx;
}


int SourceNode::x_index() const {
  pin();
  return local_->idx.x();
}


int SourceNode::y_index() const {
  pin();
  return local_->idx.y();
}


int SourceNode::z_index() const {
  pin();
  return local_->idx.z();
}


int SourceNode::level() const {
  pin();
  return local_->idx.level();
}


SourceNode SourceNode::child(size_t i) const {
  pin();
  return SourceNode{local_->child[i]};
}


SourceNode SourceNode::parent() const {
  pin();
  return SourceNode{local_->parent};
}


ExpansionRef SourceNode::expansion() const {
  pin();
  return ExpansionRef{local_->expansion};
}


SourceRef SourceNode::parts() const {
  pin();
  return local_->sources;
}


int SourceNode::n_parts() const {
  pin();
  return local_->sources.n();
}


Point SourceNode::low() const {
  pin();
  return local_->root_geo.low_from_index(local_->idx);
}


Point SourceNode::high() const {
  pin();
  return local_->root_geo.high_from_index(local_->idx);
}


Point SourceNode::center() const {
  pin();
  return local_->root_geo.center_from_index(local_->idx);
}


double SourceNode::size() const {
  pin();
  return local_->root_geo.size_from_level(local_->idx.level());
}


void SourceNode::set_expansion(std::unique_ptr<Expansion> expand) {
  ExpansionRef globexp = globalize_expansion(std::move(expand), HPX_HERE);
  pin();
  local_->expansion = globexp;
}


void SourceNode::pin() const {
  if (!local_ && data_ != HPX_NULL) {
    assert(hpx_gas_try_pin(data_, (void **)&local_));
  }
}

void SourceNode::unpin() const {
  if (local_ && data_ != HPX_NULL) {
    hpx_gas_unpin(data_);
    local_ = nullptr;
  }
}


/////////////////////////////////////////////////////////////////////
// TargetNode implementation
/////////////////////////////////////////////////////////////////////


TargetNode::TargetNode(DomainGeometry g, Index idx, hpx_addr_t method,
                       TargetNode *parent) {
  data_ = hpx_gas_alloc_local(1, sizeof(TargetNodeData), 0);
  if (data_ == HPX_NULL) {
    return;
  }
  local_ = nullptr;

  pin();
  local_->root_geo = g;
  local_->idx = idx;
  if (parent) {
    local_->parent = parent->data();
  } else {
    local_->parent = HPX_NULL;
  }
  for (int i = 0; i < 8; ++i) {
    local_->child[i] = HPX_NULL;
  }
  local_->expansion = ExpansionRef{0, HPX_NULL};
  local_->method = MethodRef{method};
  local_->targets = TargetRef{};
}


TargetNode::~TargetNode() {
  unpin();
}


void TargetNode::destroy() {
  unpin();
  hpx_call_sync(data_, target_node_delete_action, nullptr, 0, &data_);
}


void TargetNode::partition(hpx_addr_t parts, int n_parts, int limit,
                           ExpansionRef expand, int which_child,
                           bool same_sources_and_targets,
                           std::vector<SourceNode> consider) {
  hpx_addr_t done = hpx_lco_future_new(0);
  assert(done != HPX_NULL);

  size_t parms_size = target_node_partition_params_size(consider.size());
  TargetNodePartitionParams *parms =
      target_node_partition_params_alloc(consider.size());
  parms->done = done;
  parms->same_sources_and_targets = same_sources_and_targets;
  parms->parts = parts;
  parms->n_parts = n_parts;
  parms->n_parts_total = n_parts;
  parms->limit = limit;
  parms->expansion = expansion();
  parms->which_child = which_child;
  parms->n_consider = consider.size();
  for (size_t i = 0; i < consider.size(); ++i) {
    parms->consider[i] = consider[i].data();
  }

  hpx_call(data_, target_node_partition_action, HPX_NULL, parms, parms_size);
  free(parms);

  hpx_lco_wait(done);
  hpx_lco_delete_sync(done);
}


bool TargetNode::is_leaf() const {
  pin();

  bool retval{true};
  for (int i = 0; i < 8; ++i) {
    if (local_->child[i] != HPX_NULL) {
      retval = false;
      break;
    }
  }

  return retval;
}


DomainGeometry TargetNode::root_geo() const {
  pin();
  return local_->root_geo;
}


Index TargetNode::index() const {
  pin();
  return local_->idx;
}


int TargetNode::x_index() const {
  pin();
  return local_->idx.x();
}


int TargetNode::y_index() const {
  pin();
  return local_->idx.y();
}


int TargetNode::z_index() const {
  pin();
  return local_->idx.z();
}


int TargetNode::level() const {
  pin();
  return local_->idx.level();
}


TargetNode TargetNode::child(size_t i) const {
  pin();
  return TargetNode{local_->child[i]};
}


TargetNode TargetNode::parent() const {
  pin();
  return TargetNode{local_->parent};
}


ExpansionRef TargetNode::expansion() const {
  pin();
  return local_->expansion;
}


TargetRef TargetNode::parts() const {
  pin();
  return local_->targets;
}


int TargetNode::n_parts() const {
  pin();
  return local_->targets.n();
}


Point TargetNode::low() const {
  pin();
  return local_->root_geo.low_from_index(local_->idx);
}


Point TargetNode::high() const {
  pin();
  return local_->root_geo.high_from_index(local_->idx);
}


Point TargetNode::center() const {
  pin();
  return local_->root_geo.center_from_index(local_->idx);
}


double TargetNode::size() const {
  pin();
  return local_->root_geo.size_from_level(local_->idx.level());
}


void TargetNode::set_expansion(std::unique_ptr<Expansion> expand) {
  ExpansionRef globexp = globalize_expansion(std::move(expand), data_);
  pin();
  local_->expansion = globexp;
}


void TargetNode::collect_results(hpx_addr_t user_array, size_t phi_offset) {
  hpx_call_sync(data_, target_node_collect_results_action, nullptr, 0,
           &user_array, &phi_offset);
}


void TargetNode::pin() const {
  if (!local_ && data_ != HPX_NULL) {
    assert(hpx_gas_try_pin(data_, (void **)&local_));
  }
}


void TargetNode::unpin() const {
  if (local_ && data_ != HPX_NULL) {
    hpx_gas_unpin(data_);
    local_ = nullptr;
  }
}


} // namespace dashmm
