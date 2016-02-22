// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license. See the LICENSE file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================


#ifndef __DASHMM_TARGET_NODE_H__
#define __DASHMM_TARGET_NODE_H__


/// \file include/targetnode.h
/// \brief Target tree


#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>

#include <algorithm>
#include <vector>

#include <hpx/hpx.h>

#include "include/domaingeometry.h"
#include "include/expansionlco.h"
#include "include/index.h"
#include "include/point.h"
#include "include/targetlco.h"
#include "include/targetref.h"
#include "include/types.h"


namespace dashmm {



/// Forward declaration of Evaluator so that we can become friends
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class> class Method>
class Evaluator;


/// A node of the target tree.
///
/// This object is a reference object, storing only the HPX address of the
/// data for this node. As needed, the methods in this object will pin the
/// global data to produce a locally accessible version of the data to return
/// results in queries. As it is a reference, it can be passed by value
/// wherever it is needed.
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class> class Method>
class TargetNode {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Method<Source, Target, Expansion>;

  using targetref_t = TargetRef<Target>;
  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, Method>;
  using sourcenode_t = SourceNode<Source, Target, Expansion, Method>;
  using targetnode_t = TargetNode<Source, Target, Expansion, Method>;
  using targetlco_t = TargetLCO<Source, Target, Expansion, Method>;

  /// Construct from basic node data.
  ///
  /// This will create a node in the global address space with the specified
  /// data. Note that the expansion for the node will not be created; the
  /// expansion will have to be created in either inherit() or process().
  TargetNode(DomainGeometry g, Index idx, method_t method,
             targetnode_t *parent) {
    data_ = hpx_gas_alloc_local(1, sizeof(Data), 0);
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
    local_->expansion = expansionlco_t{HPX_NULL, -1};
    local_->method = method;
    local_->targets = targetlco_t{ };
  }

  /// Copy constructor
  ///
  /// This is needed so that there is correct reference counting for the
  /// pinned data locally.
  TargetNode(const targetnode_t &other)
      : local_{nullptr}, data_{other.data()} { }

  /// Construct from a global address
  explicit TargetNode(hpx_addr_t data)
      : local_{nullptr}, data_{data} { }

  /// Destructor
  ///
  /// This is required as the object may have pinned a global address, so we
  /// need to allow for the opportunity to unpin the data.
  ~TargetNode() {
    unpin();
  }

  /// Destroy the referred data
  ///
  /// To delete the node, one must call destroy(); the destruction of the
  /// TargetNode object only destroys the reference to the global data, not
  /// the global data itself.
  void destroy() {
    unpin();
    hpx_call_sync(data_, node_delete_, nullptr, 0, &data_);
  }

  /// Return the address of the referred data.
  hpx_addr_t data() const {return data_;}

  /// Partition the node.
  ///
  /// In addition to partitioning the domain based on the input target
  /// locations, this also schedules the work based on the method for the
  /// various expansions, and the target points. In particular, inherit()
  /// and process() from the method associated with the node will be called.
  ///
  /// This call is synchronous. All contributions will have been made to all
  /// targets when this returns to the calling function.
  ///
  /// Note that for internal reasons, the target points are passed to this
  /// routine as a global address, rather than as the data itself.
  ///
  /// \param targets - a reference to the targets that are in this node
  /// \param limit - the refinement limit for the partition
  /// \param n_digits - the accuracy of the expansion
  /// \param which_child - which child of this node's parent was this node
  /// \param same_sources_and_targets - are we in the case where the sources
  ///                 and targets are the same?
  /// \param consider - the list of source nodes in the consider list
  void partition(targetref_t targets, int limit,
                 int n_digits, int which_child,
                 bool same_sources_and_targets,
                 std::vector<sourcenode_t> consider) {
    hpx_addr_t done = hpx_lco_future_new(0);
    assert(done != HPX_NULL);

    size_t parms_size = partition_params_size(consider.size());
    PartitionParams *parms = partition_params_alloc(consider.size());
    parms->done = done;
    parms->same_sources_and_targets = same_sources_and_targets;
    parms->targets = targets;
    parms->limit = limit;
    parms->n_digits = n_digits;
    parms->which_child = which_child;
    parms->n_consider = consider.size();
    for (size_t i = 0; i < consider.size(); ++i) {
      parms->consider[i] = consider[i].data();
    }

    hpx_call(data_, partition_, HPX_NULL, parms, parms_size);
    delete [] parms;

    hpx_lco_wait(done);
    hpx_lco_delete_sync(done);
  }

  /// Is the node a leaf?
  bool is_leaf() const {
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

  /// Is the reference referring to data?
  bool is_valid() const {return data_ != HPX_NULL;}

  /// Return the root geometry of this node's tree
  DomainGeometry root_geo() const {
    pin();
    return local_->root_geo;
  }

  /// Return the index of this node
  Index index() const {
    pin();
    return local_->idx;
  }

  /// Return the x index of this node
  int x_index() const {
    pin();
    return local_->idx.x();
  }

  /// Return the y index of this node
  int y_index() const {
    pin();
    return local_->idx.y();
  }

  /// Return the z index of this node
  int z_index() const {
    pin();
    return local_->idx.z();
  }

  /// Return the level of this node
  int level() const {
    pin();
    return local_->idx.level();
  }

  /// Return a reference to a child of this node
  targetnode_t child(size_t i) const {
    pin();
    return targetnode_t{local_->child[i]};
  }

  /// Return a reference to the parent of this node
  targetnode_t parent() const {
    pin();
    return targetnode_t{local_->parent};
  }

  /// Return a reference to this node's expansion
  expansionlco_t expansion() const {
    pin();
    return local_->expansion;
  }

  /// Return a reference to this node's target points.
  ///
  /// If this is an internal node, the returned reference will be invalid.
  /// If this is a leaf node, the returned reference will be valid.
  targetlco_t parts() const {
    pin();
    return local_->targets;
  }

  /// Return the number of targets in this node
  ///
  /// If this is a leaf node, this will return the number of targets contained
  /// in this node. Otherwise, this will return 0.
  int n_parts() const {
    pin();
    return local_->targets.n();
  }

  /// Return the low corner of the node
  Point low() const {
    pin();
    return local_->root_geo.low_from_index(local_->idx);
  }

  /// Return the high corner of the node
  Point high() const {
    pin();
    return local_->root_geo.high_from_index(local_->idx);
  }

  /// Return the center of the node
  Point center() const {
    pin();
    return local_->root_geo.center_from_index(local_->idx);
  }

  /// Returnd the side length of this node.
  double size() const {
    pin();
    return local_->root_geo.size_from_level(local_->idx.level());
  }

  /// Set the expansion of this node.
  ///
  /// This will set this node's expansion to be a reference to the globalized
  /// form of the given expansion.
  void set_expansion(std::unique_ptr<expansion_t> expand) {
    pin();
    local_->expansion = expansionlco_t{std::move(expand), HPX_HERE};
  }

  // TODO likely this is no longer needed
  /// Collect results into user's array
  void collect_results(hpx_addr_t user_array, size_t phi_offset);

 private:
  friend class Evaluator<Source, Target, Expansion, Method>;

  /// The data for target nodes.
  ///
  /// The data stored for TargetNode objects. This will be saved in a block in
  /// the GAS.
  struct Data {
    /// The geometry of the root for the tree of which this node is a part.
    DomainGeometry root_geo;
    /// The index giving which subdivision of the root this node is.
    Index idx;
    /// The global address of the parent of this node.
    hpx_addr_t parent;
    /// The global addresses of the children of this node.
    hpx_addr_t child[8];
    /// The global address of the expansion object for this node.
    expansionlco_t expansion;
    /// The global address of the method object for this node.
    method_t method;
    /// A reference to the targets for this node.
    targetlco_t targets;
  };

  struct PartitionParams {
    hpx_addr_t done;
    // This should be expanded upon. This will be true only if the exact same
    // array is sent into evaluate as the source and target points. What this
    // means is that, if true, we can skip actual partitioning and instead
    // just find the partition point. This will then allow for a faster
    // partitioning for the target tree.
    bool same_sources_and_targets;
    targetref_t targets;
    int limit;
    int n_digits;
    int which_child;
    int n_consider;
    hpx_addr_t consider[];
  };

  void pin() const {
    if (!local_ && data_ != HPX_NULL) {
      assert(hpx_gas_try_pin(data_, (void **)&local_));
    }
  }

  void unpin() const {
    if (local_ && data_ != HPX_NULL) {
      hpx_gas_unpin(data_);
      local_ = nullptr;
    }
  }

  ///////////////////////////////////////////////////////////////////
  // Action implementations
  ///////////////////////////////////////////////////////////////////

  static int self_delete_handler(hpx_addr_t gate) {
    if (gate != HPX_NULL) {
      hpx_lco_delete_sync(gate);
    }
    hpx_addr_t target = hpx_thread_current_target();
    hpx_gas_free_sync(target);

    return HPX_SUCCESS;
  }

  static int node_delete_handler(hpx_addr_t data) {
    Data *local{nullptr};
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
          hpx_call(local->child[i], node_delete_, done, &local->child[i]);
        }
      }
    }

    hpx_gas_unpin(data);

    if (done == HPX_NULL) {
      return hpx_call_cc(data, self_delete_, &done);
    } else {
      return hpx_call_when_cc(done, data, self_delete_, &done);
    }
  }

  static size_t partition_params_size(int n_consider) {
    return sizeof(PartitionParams) + n_consider * sizeof(hpx_addr_t);
  }

  static PartitionParams *partition_params_alloc(int n_consider) {
    PartitionParams *retval = reinterpret_cast<PartitionParams *>(
      new char [partition_params_size(n_consider)]);
    if (retval) {
      retval->n_consider = n_consider;
    }
    return retval;
  }

  static int partition_handler(Data *node, PartitionParams *parms,
                               size_t bytes) {
    targetnode_t curr{hpx_thread_current_target()};
    // TODO improve this somehow.
    std::vector<sourcenode_t> consider{ };
    for (int i = 0; i < parms->n_consider; ++i) {
      consider.push_back(sourcenode_t{parms->consider[i]});
    }

    bool refine = false;
    if (parms->targets.n() > parms->limit) {
      refine = node->method.refine_test(parms->same_sources_and_targets, curr,
                                        consider);
    }

    if (!refine) {
      node->targets = targetlco_t(parms->targets);

      hpx_call_when(node->targets.lco(), parms->done, hpx_lco_set_action,
                    HPX_NULL, nullptr, 0);
    }

    node->method.inherit(curr, parms->n_digits, parms->which_child);
    node->method.process(curr, consider, !refine);

    if (refine) {
      // partition
      Target *T{nullptr};
      assert(hpx_gas_try_pin(parms->targets.data(), (void **)&T));
      Target *splits[9] { };
      splits[0] = T;
      splits[8] = &T[parms->targets.n()];

      Point cen{node->root_geo.center_from_index(node->idx)};
      double z_center = cen.z();
      auto z_comp = [&z_center](Target &a) {
                      return a.position.z() < z_center;
                    };
      double y_center = cen.y();
      auto y_comp = [&y_center](Target &a) {
                      return a.position.y() < y_center;
                    };
      double x_center = cen.x();
      auto x_comp = [&x_center](Target &a) {
                      return a.position.x() < x_center;
                    };

      if (parms->same_sources_and_targets) {
        // NOTE: The passed in value of same_sources_and_targets had better
        // be correct, or there will be trouble.
        splits[4] = std::partition_point(splits[0], splits[8], z_comp);

        splits[2] = std::partition_point(splits[0], splits[4], y_comp);
        splits[6] = std::partition_point(splits[4], splits[8], y_comp);

        splits[1] = std::partition_point(splits[0], splits[2], x_comp);
        splits[3] = std::partition_point(splits[2], splits[4], x_comp);
        splits[5] = std::partition_point(splits[4], splits[6], x_comp);
        splits[7] = std::partition_point(splits[6], splits[8], x_comp);
      } else {
        splits[4] = std::partition(splits[0], splits[8], z_comp);

        splits[2] = std::partition(splits[0], splits[4], y_comp);
        splits[6] = std::partition(splits[4], splits[8], y_comp);

        splits[1] = std::partition(splits[0], splits[2], x_comp);
        splits[3] = std::partition(splits[2], splits[4], x_comp);
        splits[5] = std::partition(splits[4], splits[6], x_comp);
        splits[7] = std::partition(splits[6], splits[8], x_comp);
      }

      hpx_gas_unpin(parms->targets.data());

      targetref_t cparts[8] {};
      int n_children{0};
      {
        int n_offset{0};
        for (int i = 0; i < 8; ++i) {
          int n_per_child = splits[i + 1] - splits[i];
          if (n_per_child) {
            ++n_children;
            cparts[i] = parms->targets.slice(n_offset, n_per_child);
          }
          n_offset += n_per_child;
        }
      }

      hpx_addr_t cdone = hpx_lco_and_new(n_children);
      assert(cdone != HPX_NULL);

      // set up the arguments to the partition actions; the constant parts
      PartitionParams *args = partition_params_alloc(consider.size());
      size_t argssize = partition_params_size(consider.size());
      args->done = cdone;
      args->same_sources_and_targets = parms->same_sources_and_targets;
      args->limit = parms->limit;
      args->n_digits = parms->n_digits;
      args->n_consider = consider.size();
      for (size_t i = 0; i < consider.size(); ++i) {
        args->consider[i] = consider[i].data();;
      }

      for (int i = 0; i < 8; ++i) {
        if (!cparts[i].valid()) {
          node->child[i] = HPX_NULL;
          continue;
        }

        Index cidx{node->idx.child(i)};
        targetnode_t thisnode{hpx_thread_current_target()};
        targetnode_t kid{node->root_geo, cidx, node->method, &thisnode};
        node->child[i] = kid.data();

        args->targets = cparts[i];
        args->which_child = i;

        hpx_call(kid.data(), partition_, HPX_NULL, args, argssize);
      }
      delete [] args;

      hpx_call_when(cdone, cdone, hpx_lco_delete_action, parms->done,
                    nullptr, 0);
    }

    // At this point, all work on the current expansion will have been scheduled,
    // so we mark the LCO as such.
    node->expansion.finalize();
    // Also, all the work on the targets for this node will have been scheduled as
    // well.
    if (!refine) {
      node->targets.finalize();
    }

    return HPX_SUCCESS;
  }

  // The actual data for the object
  mutable Data *local_;
  hpx_addr_t data_;

  // Actions for this class
  static hpx_action_t node_delete_;
  static hpx_action_t self_delete_;
  static hpx_action_t partition_;
};


template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t TargetNode<S, T, E, M>::node_delete_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t TargetNode<S, T, E, M>::self_delete_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t TargetNode<S, T, E, M>::partition_ = HPX_ACTION_NULL;


} // namespace dashmm


#endif // __DASHMM_TARGET_NODE_H__
