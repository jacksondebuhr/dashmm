#include "include/node.h"

#include <cassert.h>
#include <cstddef.h>
#include <cstdint.h>

#include "hpx/hpx.h"

#include "include/domaingeometry.h"
#include "include/index.h"
#include "include/method.h"


namespace dashmm {


struct NodeData {
  DomainGeometry root_geo;
  Index idx;
  hpx_addr_t parent;
  hpx_addr_t child[8];

  //TODO: Is this just a future that stores the ObjectBase?
  // then once it is ready, the future is set, and those reading it
  // can go ahead and get the information? I think so...
  hpx_addr_t expansion;
  hpx_addr_t method;

  //For Source Nodes, these refer to the source points
  // For target nodes, these refer to the target points
  hpx_addr_t parts;
  int n_parts;
  int n_parts_total;
};


/////////////////////////////////////////////////////////////////////
// Some shared utilities
/////////////////////////////////////////////////////////////////////


NodeData *node_data_pin(hpx_addr_t data) {
  NodeData *local{nullptr};
  assert(hpx_gas_try_pin(data, (void **)&local));
  return local;
}


void node_data_unpin(hpx_addr_t data) {
  hpx_gas_unpin(data);
}


static DomainGeometry node_get_root_geo(hpx_addr_t data) {
  NodeData *local = node_data_pin(data);
  DomainGeometry retval{local->root_geo};
  node_data_unpin(data);
  return retval;
}


static Index node_get_index(hpx_addr_t data) {
  NodeData *local = node_data_pin(data);
  Index retval{local->idx};
  node_data_unpin(data);
  return retval;
}


static int node_get_x_index(hpx_addr_t data) {
  NodeData *local = node_data_pin(data);
  Index retval{local->idx};
  node_data_unpin(data);
  return retval.x();
}


static int node_get_y_index(hpx_addr_t data) {
  NodeData *local = node_data_pin(data);
  Index retval{local->idx};
  node_data_unpin(data);
  return retval.y();
}


static int node_get_z_index(hpx_addr_t data) {
  NodeData *local = node_data_pin(data);
  Index retval{local->idx};
  node_data_unpin(data);
  return retval.z();
}


static int node_get_level(hpx_addr_t data) {
  NodeData *local = node_data_pin(data);
  Index retval{local->idx};
  node_data_unpin(data);
  return retval.level();
}


static hpx_addr_t node_get_child(hpx_addr_t data, size_t i) {
  NodeData *local = node_data_pin(data);
  hpx_addr_t retval{local->child[i]};
  node_data_unpin(data);
  return retval;
}


static hpx_addr_t node_get_parent(hpx_addr_t data) {
  NodeData *local = node_data_pin(data);
  hpx_addr_t retval{local->parent};
  node_data_unpin(data);
  return retval;
}


static hpx_addr_t node_get_expansion(hpx_addr_t data) {
  NodeData *local = node_data_pin(data);
  hpx_addr_t retval{local->expansion};
  node_data_unpin(data);
  return retval;
}


static hpx_addr_t node_get_parts(hpx_addr_t data) {
  NodeData *local = node_data_pin(data);
  hpx_addr_t retval{local->parts};
  node_data_unpin(data);
  return retval;
}


static int node_get_n_parts(hpx_addr_t data) {
  NodeData *local = node_data_pin(data);
  int retval{local->n_parts};
  node_data_unpin(data);
  return retval;
}


static int node_get_n_parts_total(hpx_addr_t data) {
  NodeData *local = node_data_pin(data);
  int retval{local->n_parts_total};
  node_data_unpin(data);
  return retval;
}


static Point node_get_low(hpx_addr_t data) {
  NodeData *local = node_data_pin(data);
  Point retval = local->geo.low_from_index(local->idx.x(), local->idx.y,
                                           local->idx.z(), local->idx.level());
  node_data_unpin(data);
  return retval;
}


static Point node_get_high(hpx_addr_t data) {
  NodeData *local = node_data_pin(data);
  Point retval = local->geo.high_from_index(local->idx.x(), local->idx.y,
                                           local->idx.z(), local->idx.level());
  node_data_unpin(data);
  return retval;
}


static Point node_get_center(hpx_addr_t data) {
  NodeData *local = node_data_pin(data);
  Point retval = local->geo.center_from_index(local->idx.x(), local->idx.y,
                                           local->idx.z(), local->idx.level());
  node_data_unpin(data);
  return retval;
}


static Point node_get_size(hpx_addr_t data) {
  NodeData *local = node_data_pin(data);
  Point retval = local->geo.size_from_index(local->idx.level());
  node_data_unpin(data);
  return retval;
}


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


int node_delete_handler(hpx_addr_t data) {
  NodeData *local{nullptr};
  assert(hpx_gas_try_pin(data, (void **)&local));

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
        hpx_call(local->child[i], node_delete_action, done,
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
HPX_ACTION(HPX_DEFAULT, 0, node_delete_action, node_delete_handler,
           HPX_ADDR_T);


int source_node_child_generation_done_handler(NodeData *node,
                                              hpx_addr_t gendone,
                                              hpx_addr_t expand) {
  ExpansionRef expref{expand};
  MethodRef method{node->method};
  SourceNode snode{hpx_thread_current_target()};

  method.aggregate(snode, expref);

  hpx_lco_delete_sync(gendone);
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           source_node_child_generation_done_action,
           source_node_child_generation_done_handler,
           HPX_POINTER, HPX_ADDR, HPX_ADDR);


int source_node_partition_handler(NodeData *node, hpx_addr_t partdone,
                                  hpx_addr_t gendone, hpx_addr_t parts,
                                  int n_parts, int n_parts_total, int limit,
                                  hpx_addr_t expand);
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           source_node_partition_action, source_node_partition_handler,
           HPX_POINTER, HPX_ADDR, HPX_ADDR, HPX_ADDR, HPX_INT, HPX_INT,
           HPX_INT, HPX_ADDR);

int source_node_partition_handler(NodeData *node, hpx_addr_t partdone,
                                  hpx_addr_t gendone, hpx_addr_t parts,
                                  int n_parts, int n_parts_total, int limit,
                                  hpx_addr_t expand) {
  assert(n_parts > 0);
  //NOTE: The following line only makes sense in SMP; In distributed, having
  // the address of the start of parts will not work, as those are temporary
  // global addresses. I think.
  node->parts = parts;
  node->n_parts = n_parts;
  node->n_parts_total = n_parts_total;

  if (n_parts <= limit) {
    MethodRef method{node->method};
    SourceNode curr{hpx_thread_current_target()};
    ExpansionRef expref{expand};
    method.generate(curr, expref);

    hpx_lco_set(alldone, 0, nullptr, HPX_NULL, HPX_NULL);
    return HPX_SUCCESS;
  }

  //NOTE: We make use of the SMP only form here. This will not work in
  // distributed.

  //partition sources
  Source *P{nullptr};
  assert(hpx_gas_try_pin(parts, (void **)&P));
  Source *splits[9]{};
  splits[0] = P;
  splits[8] = &P[n_parts];

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
  hpx_gas_unpin(parts);

  //copy sections of particles into their own allocations
  // This will be needed in distributed mode.
  hpx_addr_t cparts[8];
  int n_per_child[8]{0, 0, 0, 0, 0, 0, 0, 0};
  int n_children{0};
  int n_offset{0};
  for (int i = 0; i < 8; ++i) {
    n_per_child[i] = splits[i + 1] - splits[i];
    if (n_per_child[i]) {
      ++n_children;
      cparts[i] = hpx_addr_add(parts, sizeof(Source) * n_offset,
                                      sizeof(Source) * n_parts_total);
    } else {
      cparts[i] = HPX_NULL;
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

    hpx_call(kid.data(), source_node_partition_action, HPX_NULL, &childpartdone,
             &childgendone, &cparts[i], &n_per_child[i], &n_parts_total, &limit,
             &expand);
  }

  hpx_call_when(childpartdone, childpartdone, hpx_lco_delete_action,
                partdone, nullptr, 0);
  hpx_call_when(childgendone, hpx_thread_current_target(),
                source_node_child_generation_done_action, gendone,
                &childgendone, &expand);

  return HPX_NULL;
}


struct TargetNodePartitionParams {
  hpx_addr_t done;
  bool same_sources_and_targets;
  hpx_addr_t parts;
  int n_parts;
  int n_parts_total;
  int limit;
  hpx_addr_t expansion;
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
  node->parts = parms->parts;
  node->n_parts = parms->n_parts;
  node->n_parts_total = parms->n_parts_total;

  bool refine = false;
  MethodRef method{node->method};
  TargetNode curr{hpx_thread_current_target()};
  std::vector<SourceNode> consider{};
  for (size_t i = 0; i < parms->n_consider; ++i) {
    consider.push_back(SourceNode{parms->consider[i]});
  }

  if (parms->n_parts > parms->limit) {
    refine = method.refine_test(parms->same_sources_and_targets, curr,
                                consider);
  }

  ExpansionRef expand{parms->expand};
  method.inherit(curr, expand, parms->which_child);
  //If we are not going to refine, then curr is a leaf
  method.process(curr, consider, !refine);

  if (refine) {
    //partition
    Target *P{nullptr};
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
      TargetNode kid{node->root_geo, cidx.x(), cidx.y(), cidx.z(), cidx.level(),
                     node->method, hpx_thread_current_target()};
      node->child[i] = kid.data();

      args->parts = cparts[i];
      args->n_parts = n_per_child[i];
      args->which_child = i;

      hpx_call(kid.data(), target_node_partition_action, HPX_NULL, args,
               argssize);
    }

    hpx_call_when(cdone, cdone, hpx_lco_delete_action, parms->done, nullptr, 0);
  } else {  // no refinement needed; so just set the input done LCO
    hpx_lco_and_set(parms->done, HPX_NULL);
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

  NodeData *local = node_data_pin(data_);
  local->root_geo = g;
  local->idx = Index{ix, iy, iz, level};
  local->parent = parent->data();
  for (int i = 0; i < 8; ++i) {
    local->child[0] = HPX_NULL;
  }
  local->expansion = HPX_NULL;
  local->method = method;
  local->parts = HPX_NULL;
  local->n_parts = 0;
  node_data_unpin(data_);
}


SourceNode::~SourceNode() {
  hpx_call_sync(data_, node_delete_action, nullptr, 0, &data_);
}


hpx_addr_t SourceNode::partition(hpx_addr_t parts, int n_parts, int limit,
                                 hpx_addr_t expand) {
  hpx_addr_t retval = hpx_lco_future_new(0);
  assert(retval != HPX_NULL);

  hpx_addr_t empty{HPX_NULL};
  hpx_call(data_, source_node_partition_action, HPX_NULL,
           &retval, &empty, &parts, &n_parts, &limit, &expand);

  return retval;
}


bool SourceNode::is_leaf() const {
  NodeData *local = node_data_pin(data_);
  bool retval{true};
  for (int i = 0; i < 8; ++i) {
    if (child[i] != HPX_NULL) {
      retval = false;
      break;
    }
  }
  node_data_unpin(data_);
  return retval;
}


DomainGeometry SourceNode::root_geo() const {
  return node_get_root_geo(data_);
}


Index SourceNode::index() const {
  return node_get_index(data_);
}


int SourceNode::x_index() const {
  return node_get_x_index(data_);
}


int SourceNode::y_index() const {
  return node_get_y_index(data_);
}


int SourceNode::z_index() const {
  return node_get_z_index(data_);
}


int SourceNode::level() const {
  return node_get_level(data_);
}


SourceNode *SourceNode::child(size_t i) const {
  return SourceNode{node_get_child(data_, i)};
}


SourceNode *SourceNode::parent() const {
  return SourceNode{node_get_parent(data_)};
}


ExpansionRef SourceNode::expansion() const {
  return ExpansionRef{node_get_expansion(data_)};
}


SourceRef SourceNode::parts() const {
  NodeData *local = node_data_pin(data);
  SourceRef retval{local->parts, local->n_parts, local->n_parts_total};
  node_data_unpin(data);
  return retval;
}


int SourceNode::n_parts() const {
  return node_get_n_parts(data_);
}


int SourceNode::n_parts_total() const {
  return node_get_n_parts_total(data_);
}


Point SourceNode::low() const {
  return node_get_low(data_);
}


Point SourceNode::high() const {
  return node_get_high(data_);
}


Point SourceNode::center() const {
  return node_get_center(data_);
}


double SourceNode::size() const {
  return node_get_size(data_);
}


void SourceNode::set_expansion(ExpansionRef expand) {
  NodeData *local = node_data_pin(data_);
  local->expansion = expand.data();
  node_data_unpin(data_);
}


/////////////////////////////////////////////////////////////////////
// TargetNode implementation
/////////////////////////////////////////////////////////////////////


TargetNode::TargetNode(DomainGeometry g, int ix, int iy, int iz, int level,
                       hpx_addr_t method, TargetNode *parent) {
  data_ = hpx_gas_alloc_local(sizeof(NodeData), 0);
  if (data_ == HPX_NULL) {
    return;
  }

  NodeData *local = node_data_pin(data_);
  local->root_geo = g;
  local->idx = Index{ix, iy, iz, level};
  local->parent = parent->data();
  for (int i = 0; i < 8; ++i) {
    local->child[0] = HPX_NULL;
  }
  local->expansion = HPX_NULL;
  local->method = method;
  local->parts = HPX_NULL;
  local->n_parts = 0;
  node_data_unpin(data_);
}


TargetNode::~TargetNode() {
  hpx_call_sync(data_, node_delete_action, nullptr, 0, &data_);
}


void TargetNode::partition(hpx_addr_t parts, int n_parts, int limit,
                           hpx_addr_t expand, int which_child,
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
  NodeData *local = node_data_pin(data_);
  bool retval{true};
  for (int i = 0; i < 8; ++i) {
    if (child[i] != HPX_NULL) {
      retval = false;
      break;
    }
  }
  node_data_unpin(data_);
  return retval;
}


DomainGeometry TargetNode::root_geo() const {
  return node_get_root_geo(data_);
}


Index TargetNode::index() const {
  return node_get_index(data_);
}


int TargetNode::x_index() const {
  return node_get_x_index(data_);
}


int TargetNode::y_index() const {
  return node_get_y_index(data_);
}


int TargetNode::z_index() const {
  return node_get_z_index(data_);
}


int TargetNode::level() const {
  return node_get_level(data_);
}


TargetNode *TargetNode::child(size_t i) const {
  return TargetNode{node_get_child(data_, i)};
}


TargetNode *TargetNode::parent() const {
  return TargetNode{node_get_parent(data_, i)};
}


ExpansionRef TargetNode::expansion() const {
  return ExpansionRef{node_get_expansion(data_)};
}


TargetRef TargetNode::parts() const {
  NodeData *local = node_data_pin(data);
  TargetRef retval{local->parts, local->n_parts, local->n_parts_total};
  node_data_unpin(data);
  return retval;
}


int TargetNode::n_parts() const {
  return node_get_n_parts(data_);
}


int TargetNode::n_parts_total() const {
  return node_get_n_parts_total(data_);
}


Point TargetNode::low() const {
  return node_get_low(data_);
}


Point TargetNode::high() const {
  return node_get_high(data_);
}


Point TargetNode::center() const {
  return node_get_center(data_);
}


double TargetNode::size() const {
  return node_get_size(data_);
}


void TargetNode::set_expansion(ExpansionRef expand) {
  NodeData *local = node_data_pin(data_);
  local->expansion = expand.data();
  node_data_unpin(data_);
}


} // namespace dashmm
