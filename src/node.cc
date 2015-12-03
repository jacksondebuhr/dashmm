#include "include/node.h"

#include <cassert.h>
#include <cstddef.h>
#include <cstdint.h>

#include "hpx/hpx.h"

#include "include/domaingeometry.h"
#include "include/index.h"
#include "include/method.h"


namespace dashmm {


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
           HPX_ADDR_T);


int source_node_delete_handler(hpx_addr_t data) {
  SourceNodeData *local{nullptr};
  assert(hpx_gas_try_pin(data, (void **)&local));

  local->expansion.destroy();
  //TODO: If we really are making copies of the method always, go ahead and
  // destoy it here. Otherwise, there is a single instance of it somewhere.
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
HPX_ACTION(HPX_DEFAULT, 0,
           source_node_delete_action, source_node_delete_handler,
           HPX_ADDR_T);


int target_node_delete_handler(hpx_addr_t data) {
  TargetNodeData *local{nullptr};
  assert(hpx_gas_try_pin(data, (void **)&local));

  local->expansion.destroy();
  //TODO: If we really are making copies of the method always, go ahead and
  // destoy it here. Otherwise, there is a single instance of it somewhere.
  local->targets.destroy();
  hpx_lco_delete_sync(local->part_done);

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
HPX_ACTION(HPX_DEFAULT, 0,
           target_node_delete_action, target_node_delete_handler,
           HPX_ADDR_T);


int source_node_child_generation_done_handler(NodeData *node,
                                              hpx_addr_t gendone,
                                              hpx_addr_t expand) {
  ExpansionRef expref{expand};
  MethodRef method{node->method};
  SourceNode snode{hpx_thread_current_target()};

  method.aggregate(snode, expref);

  //TODO: set the all scheduled input on the expansion

  hpx_lco_delete_sync(gendone);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           source_node_child_generation_done_action,
           source_node_child_generation_done_handler,
           HPX_POINTER, HPX_ADDR, HPX_ADDR);


struct SourceNodePartitionParams {
  hpx_addr_t partdone;
  hpx_addr_t gendone;
  int limit;
  int type;
  hpx_addr_t expand;
  int n_sources;
  Source sources[];
};

int source_node_partition_handler(NodeData *node,
                                  SourceNodePartitionParams *parms,
                                  size_t bytes);
HPX_ACTION(HPX_DEFAULT, HPX_PINNED | HPX_MARSHALLED,
           source_node_partition_action, source_node_partition_handler,
           HPX_POINTER, HPX_POINTER, HPX_SIZE_T);

int source_node_partition_handler(NodeData *node,
                                  SourceNodePartitionParams *parms,
                                  size_t bytes) {
  if (parms->n_sources <= parms->limit) {
    assert(node->sources.data() == HPX_NULL);
    node->sources = SourceRef(parms->sources, parms->n_sources);

    MethodRef method{node->method};
    SourceNode curr{hpx_thread_current_target()};
    ExpansionRef expref{parms->type, parms->expand};
    method.generate(curr, expref);

    //TODO: set the all scheduled input on the expansion

    hpx_lco_set(parms->partdone, 0, nullptr, HPX_NULL, HPX_NULL);
    if (parms->gendone != HPX_NULL) {
      hpx_lco_set(parms->gendone, 0, nullptr, HPX_NULL, HPX_NULL);
    }
    return HPX_SUCCESS;
  }

  //NOTE: We make use of the SMP only form here. This will not work in
  // distributed.

  //partition sources
  Source *splits[9]{};
  splits[0] = parms->sources;
  splits[8] = &parms->sources[n_parts];

  Point cen{center()};
  double z_center = cen.z();
  auto z_comp = [&z_center](Source a) {
    return a.z() < z_center;
  };
  splits[4] = std::partition(splits[0], splits[8], z_comp);

  double y_center = cen.y();
  auto y_comp = [&y_center](Source a) {
    return a.y() < y_center;
  };
  splits[2] = std::partition(splits[0], splits[4], y_comp);
  splits[6] = std::partition(splits[4], splits[8], y_comp);

  double x_center = cen.x();
  auto x_comp = [&x_center](Source a) {
    return a.x() < x_center;
  };
  splits[1] = std::partition(splits[0], splits[2], x_comp);
  splits[3] = std::partition(splits[2], splits[4], x_comp);
  splits[5] = std::partition(splits[4], splits[6], x_comp);
  splits[7] = std::partition(splits[6], splits[8], x_comp);

  //copy sections of particles into their own allocations
  // This will be needed in distributed mode.
  Source *cparts[8];
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
  hpx_addr_t childgendone = hpx_lco_and_new(n_children);
  assert(childgendone != HPX_NULL);

  for (int i = 0; i < 8; ++i) {
    if (n_per_child[i] == 0) {
      node->child[i] = HPX_NULL;
      continue;
    }

    Index cidx{node->idx.child(i)};
    SourceNode kid{node->root_geo, cidx.x(), cidx.y(), cidx.z(), cidx.level(),
                   node->method, hpx_thread_current_target()};
    node->child[i] = kid.data();

    //TODO: can probably pretty easily remove this allocation/deallocation
    // overhead.
    size_t argsz = sizeof(SourceNodePartitionParams)
                   + sizeof(Sourc) * n_per_child[i]);
    SourceNodePartitionParams *args = malloc(argsz);
    assert(args != nullptr);
    args->partdone = childpartdone;
    args-gendone = childgendone;
    args->limit = parms->limit;
    args->type = parms->type;
    args->expand = parms->expand;
    args->n_sources = n_per_child[i];
    memcpy(args->sources, cparts[i], sizeof(Source) * n_per_child[i]);
    hpx_call(kid.data(), source_node_partition_action, HPX_NULL, args, argsz);
    free(args);
  }

  hpx_call_when(childpartdone, childpartdone, hpx_lco_delete_action,
                parms->partdone, nullptr, 0);
  hpx_call_when(childgendone, hpx_thread_current_target(),
                source_node_child_generation_done_action, parms->gendone,
                &childgendone, &parms->expand);

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
      malloc(target_node_partition_params_size(n_consider));
  if (retval) {
    retval->n_consider = n_consider;
  }
  return retval;
}

int target_node_partition_handler(NodeData *node,
                                  TargetNodePartitionParams *parms,
                                  size_t bytes);
HPX_ACTION(HPX_DEFAULT, HPX_PINNED | HPX_MARSHALLED,
           target_node_partition_action, target_node_partition_handler,
           HPX_POINTER, HPX_POINTER, HPX_SIZE_T);

int target_node_partition_handler(NodeData *node,
                                  TargetNodePartitionParams *parms,
                                  size_t bytes) {
  MethodRef method{node->method};
  TargetNode curr{hpx_thread_current_target()};
  std::vector<SourceNode> consider{};
  for (size_t i = 0; i < parms->n_consider; ++i) {
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
    assert(hpx_gas_try_pin(parms->parts, (void **)T));
    Target *splits[9]{};
    splits[0] = T;
    splits[8] = &T[parms->n_parts];

    Point cen{center()};
    double z_center = cen.z();
    auto z_comp = [&z_center](Target a) {
      return a.z() < z_center;
    };
    splits[4] = std::partition(splits[0], splits[8], z_comp);

    double y_center = cen.y();
    auto y_comp = [&y_center](Target a) {
      return a.y() < y_center;
    };
    splits[2] = std::partition(splits[0], splits[4], y_comp);
    splits[6] = std::partition(splits[4], splits[8], y_comp);

    double x_center = cen.x();
    auto x_comp = [&x_center](Target a) {
      return a.x() < x_center;
    };
    splits[1] = std::partition(splits[0], splits[2], z_comp);
    splits[3] = std::partition(splits[2], splits[4], z_comp);
    splits[5] = std::partition(splits[4], splits[6], z_comp);
    splits[7] = std::partition(splits[6], splits[8], z_comp);
    hpx_gas_unpin(parms->parts);

    hpx_addr_t cparts[8];
    int n_per_child[8]{0, 0, 0, 0, 0, 0, 0, 0};
    int n_children{0};
    int n_offset{0};
    for (int i = 0; i < 8; ++i) {
      n_per_child[i] = splits[i + 1] - splits[i];
      if (n_per_child[i]) {
        ++n_children;
        cparts[i] = hpx_addr_add(parts, sizeof(Source) * n_offset,
                                        sizeof(Source) * parms->n_parts_total);
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
      TargetNode kid{node->root_geo, cidx, node->method,
                     hpx_thread_current_target()};
      node->child[i] = kid.data();

      args->parts = cparts[i];
      args->n_parts = n_per_child[i];
      args->which_child = i;

      hpx_call(kid.data(), target_node_partition_action, HPX_NULL, args,
               argssize);
    }

    hpx_call_when(cdone, node->part_done, hpx_lco_set_action, parms->done,
                  nullptr, 0);
    hpx_call_when(cdone, cdone, hpx_lco_delete_action, HPX_NULL, nullptr, 0);
  } else {  // no refinement needed; so just set the input done LCO
    Targets *targs{nullptr};
    assert(hpx_gas_try_pin(parms->parts, (void **)&targs));
    node->targets = TargetRef(targs, parms->n_parts);

    hpx_lco_and_set(parms->done, HPX_NULL);
    hpx_lco_set(node->part_done, 0, nullptr, HPX_NULL, HPX_NULL);
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


/////////////////////////////////////////////////////////////////////
// SourceNode implementation
/////////////////////////////////////////////////////////////////////


SourceNode::SourceNode(DomainGeometry g, int ix, int iy, int iz, int level,
                       hpx_addr_t method, SourceNode *parent) {
  data_ = hpx_gas_alloc_local(sizeof(NodeData), 0);
  if (data_ == HPX_NULL) {
    return;
  }

  pin();
  local_->root_geo = g;
  local_->idx = Index{ix, iy, iz, level};
  local_->parent = parent->data();
  for (int i = 0; i < 8; ++i) {
    local_->child[0] = HPX_NULL;
  }
  local_->expansion = HPX_NULL;
  local_->method = method;
  local_->sources = SourceRef{};
}


SourceNode::~SourceNode() {
  unpin();
}


void SourceNode::destroy() {
  hpx_call_sync(data_, source_node_delete_action, nullptr, 0, &data_);
}


hpx_addr_t SourceNode::partition(Source *sources, int n_sources, int limit,
                                 hpx_addr_t expand) {
  hpx_addr_t retval = hpx_lco_future_new(0);
  assert(retval != HPX_NULL);

  hpx_addr_t empty{HPX_NULL};
  hpx_call(data_, source_node_partition_action, HPX_NULL,
           &retval, &empty, &parts, &n_parts, &limit, &expand);

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
  return SourceRef{local->parts, local->n_parts, local->n_parts_total};
}


int SourceNode::n_parts() const {
  pin();
  return local_->n_parts;
}


int SourceNode::n_parts_total() const {
  pin();
  return local_->n_parts_total;
}


Point SourceNode::low() const {
  pin();
  return local->root_geo.low_from_index(local_->idx.x(), local_->idx.y(),
                                   local_->idx.z(), local_->idx.level());
}


Point SourceNode::high() const {
  pin();
  return local->root_geo.high_from_index(local_->idx.x(), local_->idx.y(),
                                   local_->idx.z(), local_->idx.level());
}


Point SourceNode::center() const {
  pin();
  return local->root_geo.center_from_index(local_->idx.x(), local_->idx.y(),
                                   local_->idx.z(), local_->idx.level());
}


double SourceNode::size() const {
  pin();
  return local->root_geo.size_from_index(local_->idx.level());
}


void SourceNode::set_expansion(std::unique_ptr<Expansion> expand) {
  //TODO: update this for new interface
  ExpansionRef globexp = globalize_expansion(expand.get(), data_);
  pin();
  local_->expansion = globexp.data();
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
  data_ = hpx_gas_alloc_local(sizeof(NodeData), 0);
  if (data_ == HPX_NULL) {
    return;
  }

  pin();
  local_->root_geo = g;
  local_->idx = idx;
  local_->parent = parent->data();
  for (int i = 0; i < 8; ++i) {
    local_->child[0] = HPX_NULL;
  }
  local_->expansion = HPX_NULL;
  local_->method = method;
  local_->targets = TargetRef{};
}


TargetNode::~TargetNode() {
  unpin();
}


void TargetNode::destroy() {
  hpx_call_sync(data_, target_node_delete_action, nullptr, 0, &data_);
}


void TargetNode::partition(hpx_addr_t parts, int n_parts, int limit,
                           ExpansionRef expand, int which_child,
                           bool same_sources_and_targets,
                           std::vector<SourceNode *> consider) {
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
    parms->consider[i] = consider[i]->data();
  }

  hpx_call(data_, target_node_partition_action, HPX_NULL, parms, parms_size);

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
  return ExpansionRef{local->expansion};
}


TargetRef TargetNode::parts() const {
  pin();
  return TargetRef{local->parts, local->n_parts, local->n_parts_total};;
}


int TargetNode::n_parts() const {
  pin();
  return local_->n_parts;
}


int TargetNode::n_parts_total() const {
  pin();
  return local_->n_parts_total;
}


Point TargetNode::low() const {
  pin();
  return local_->root_geo.low_from_index(local_->idx.x(), local_->idx.y(),
                                   local_->idx.z(), local_->idx.level());
}


Point TargetNode::high() const {
  pin();
  return local_->root_geo.high_from_index(local_->idx.x(), local_->idx.y(),
                                   local_->idx.z(), local_->idx.level());
}


Point TargetNode::center() const {
  pin();
  return local_->root_geo.center_from_index(local_->idx.x(), local_->idx.y(),
                                   local_->idx.z(), local_->idx.level());
}


double TargetNode::size() const {
  pin();
  return local_->root_geo.size_from_index(local_->idx.level());
}


void TargetNode::set_expansion(std::unique_ptr<Expansion> expand) {
  //TODO: update this for new interface
  ExpansionRef globexp = globalize_expansion(expand.get(), data_);
  pin();
  local->expansion = globexp;
}


void TargetNode::pin() {
  if (!local_ && data_ != HPX_NULL) {
    assert(hpx_gas_try_pin(data_, (void **)&local));
  }
}


void TargetNode::unpin() {
  if (local_ && data_ != HPX_NULL) {
    hpx_gas_unpin(data_);
    local_ = nullptr;
  }
}


} // namespace dashmm
