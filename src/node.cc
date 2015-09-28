#include "include/node.h"

#include <cassert.h>
#include <cstddef.h>
#include <cstdint.h>

#include "hpx/hpx.h"

#include "include/domaingeometry.h"
#include "include/index.h"


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
  hpx_addr_t first;
  hpx_addr_t last;
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


static hpx_addr_t node_get_first(hpx_addr_t data) {
  NodeData *local = node_data_pin(data);
  hpx_addr_t retval{local->first};
  node_data_unpin(data);
  return retval;
}


static hpx_addr_t node_get_last(hpx_addr_t data) {
  NodeData *local = node_data_pin(data);
  hpx_addr_t retval{local->last};
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
  local->first = HPX_NULL;
  local->last = HPX_NULL;
  node_data_unpin(data_);
}


SourceNode::~SourceNode() {
  hpx_call_sync(data_, node_delete_action, nullptr, 0, &data_);
}


//TODO
hpx_addr_t SourceNode::partition(hpx_addr_t first, hpx_addr_t last, int limit,
                                 hpx_addr_t expand) {
  //create the partition signal LCO
  //call the action, passing in that LCO
  //return LCO addx
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


hpx_addt_t SourceNode::expansion() const {
  return node_get_expansion(data_);
}


hpx_addr_t SourceNode::first() const {
  return node_get_first(data_);
}


hpx_addr_t SourceNode::last() const {
  return node_get_last(data_);
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


void SourceNode::set_expansion(hpx_addr_t expand) {
  NodeData *local = node_data_pin(data_);
  local->expansion = expand;
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
  local->first = HPX_NULL;
  local->last = HPX_NULL;
  node_data_unpin(data_);
}


TargetNode::~TargetNode() {
  hpx_call_sync(data_, node_delete_action, nullptr, 0, &data_);
}


//TODO
void TargetNode::partition(std::vector<Target>::iterator first,
                           std::vector<Target>::iterator last, int limit,
                           hpx_addr_t expand, int which_child,
                           std::vector<SourceNode *> consider) {
  //
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


hpx_addr_t TargetNode::expansion() const {
  return node_get_expansion(data_);
}


hpx_addr_t TargetNode::first() const {
  return node_get_first(data_);
}


hpx_addr_t TargetNode::last() const {
  return node_get_last(data_);
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


void TargetNode::set_expansion(Expansion *expand) {
  NodeData *local = node_data_pin(data_);
  local->expansion = expand;
  node_data_unpin(data_);
}


} // namespace dashmm
