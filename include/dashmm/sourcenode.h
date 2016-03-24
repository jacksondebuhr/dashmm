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


#ifndef __DASHMM_SOURCE_NODE_H__
#define __DASHMM_SOURCE_NODE_H__


/// \file include/dashmm/sourcenode.h
/// \brief Source tree


#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>

#include <algorithm>
#include <memory>

#include <hpx/hpx.h>

#include "dashmm/domaingeometry.h"
#include "dashmm/expansionlco.h"
#include "dashmm/index.h"
#include "dashmm/point.h"
#include "dashmm/sourceref.h"
#include "dashmm/types.h"


namespace dashmm {


/// Forward declaration of Evaluator so that we can become friends
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class> class Method>
class Evaluator;


/// A node of the source tree.
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
class SourceNode {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Method<Source, Target, Expansion>;

  using sourceref_t = SourceRef<Source>;
  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, Method>;
  using sourcenode_t = SourceNode<Source, Target, Expansion, Method>;

  /// Construct from basic node data.
  SourceNode(DomainGeometry g, Index idx,
             method_t met, sourcenode_t *parent) {
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
    local_->method = met;
    local_->sources = sourceref_t{ };
  }

  /// Construct from a global address.
  explicit SourceNode(hpx_addr_t data = HPX_NULL)
      : local_{nullptr}, data_{data} { }

  /// Copy constructor
  ///
  /// This is needed to assure that pin()/unpin() are properly reference counted
  /// when new instances of SourceNode refer to the same global data. If
  /// local_ was accidentally copied from a pinned node, the unpin() on either
  /// would invalidate the local pointer for the other.
  SourceNode(const sourcenode_t &other)
      : local_{nullptr}, data_{other.data()} { }

  /// Destructor
  ///
  /// This is required because the object may have pinned a global address,
  /// and so this address will have to be unpinned.
  ~SourceNode() {
    unpin();
  }

  /// Destroy the referred to data
  ///
  /// To delete the node, one must call destroy(); the destruction of the
  /// SourceNode object only destroys the reference to the global data, not
  /// the global data itself.
  void destroy() {
    unpin();
    hpx_call_sync(data_, node_delete_, nullptr, 0, &data_);
  }

  /// Return the address of the referred to data.
  hpx_addr_t data() const {return data_;}

  /// Partition the node.
  ///
  /// In addition to partitioning the domain based on the input particles, and
  /// input partition limit, this also schedules the work based on the method
  /// and the given expansion. In particular, the generate() and aggregate()
  /// methods will be called for this node based on the leaf/internal status
  /// of this node after partitioning.
  ///
  /// Partitioning will finish before the work of the method is finished, and
  /// for internal nodes, the difference in completion time could be large
  /// because the child node may have to wait itself for moments from its
  /// children before it can aggregate its own moments.
  ///
  /// This will return the address of an LCO that represents completion of the
  /// partitioning of this node.
  ///
  /// In practice, this is only called once in Evaluator::evaluate(). The
  /// internal calls to partition are handled directly.
  ///
  /// \param sources - the source points
  /// \param limit - the refinement limit for the partitioning
  /// \param n_digits - accuracy of the expansion
  hpx_addr_t partition(sourceref_t sources, int limit, int n_digits) {
    hpx_addr_t retval = hpx_lco_future_new(0);
    assert(retval != HPX_NULL);

    PartitionParams args{ };
    args.partdone = retval;
    args.limit = limit;
    args.n_digits = n_digits;
    args.sources = sources;
    hpx_call(data_, partition_, HPX_NULL, &args, sizeof(args));

    return retval;
  }

  /// Is the referred to node a leaf?
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

  /// Does this refer to data?
  bool is_valid() const {return data_ != HPX_NULL;}

  /// Return the geometry of the root of the tree
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

  /// Return reference to a child of this node
  sourcenode_t child(size_t i) const {
    pin();
    return sourcenode_t{local_->child[i]};
  }

  /// Return reference to the parent of this node
  sourcenode_t parent() const {
    pin();
    return sourcenode_t{local_->parent};
  }

  /// Return this node's expansion
  expansionlco_t expansion() const {
    pin();
    return local_->expansion;
  }

  /// Return the source points for this node
  ///
  /// For a leaf node, this will return a reference to the source poitns in
  /// that node. For an internal node, this will return an empty reference
  sourceref_t parts() const {
    pin();
    return local_->sources;
  }

  /// Return the number of source points for this node
  ///
  /// For a leaf node, this will return the number of source points in that
  /// leaf. For an internal node, this will return 0.
  int n_parts() const {
    pin();
    return local_->sources.n();
  }

  /// Return the low corner of this node
  Point low() const {
    pin();
    return local_->root_geo.low_from_index(local_->idx);
  }

  /// Return the high corner of this node
  Point high() const {
    pin();
    return local_->root_geo.high_from_index(local_->idx);
  }

  /// Return the center of this node
  Point center() const {
    pin();
    return local_->root_geo.center_from_index(local_->idx);
  }

  /// Return the edge size of this node
  double size() const {
    pin();
    return local_->root_geo.size_from_level(local_->idx.level());
  }

  /// Set the expansion of this node.
  ///
  /// This sets the expansion of this node to be a reference to a globalized
  /// version of the probided expansion.
  ///
  /// \param expand - the expansion to set as the expansion of this node
  void set_expansion(std::unique_ptr<expansion_t> expand) {
    pin();
    local_->expansion = expansionlco_t{std::move(expand), HPX_HERE};
  }

 private:
  friend class Evaluator<Source, Target, Expansion, Method>;

  /// The data for source nodes.
  ///
  /// The data stored for SourceNode objects. This will be saved in a block in
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
    /// A reference to the sources for this node. If this is an internal node,
    /// this will be an invalid SourceRef.
    sourceref_t sources;
  };

  /// Parameters to partition action
  struct PartitionParams {
    hpx_addr_t partdone;
    size_t limit;
    int n_digits;
    sourceref_t sources;
  };

  // The current strategy for dealing with the fact that the data lives in
  // GAS is to just pin whenever. This is fine in SMP. This will need to
  // be improved upon in the future.
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
    // NOTE: We do not destroy sources, as it is a shared resource.

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

  static int child_done_handler(Data *node, hpx_addr_t partdone,
                                int n_digits) {
    sourcenode_t snode{hpx_thread_current_target()};

    node->method.aggregate(snode, n_digits);
    node->expansion.finalize();

    hpx_lco_delete_sync(partdone);
    return HPX_SUCCESS;
  }

  static int partition_handler(Data *node, PartitionParams *parms,
                               size_t bytes) {
    // Have we reached a leaf?
    if (parms->sources.n() <= parms->limit) {
      assert(node->sources.data() == HPX_NULL);
      node->sources = parms->sources;

      sourcenode_t curr{hpx_thread_current_target()};

      // This will cause the creation of a new expansion, that will be globalized
      // and then set as the expansion for this node. It will be created with
      // initial data, and so there is no need to schedule any contributions to
      // the expansion.
      node->method.generate(curr, parms->n_digits);

      // We are done scheduling contributions to the expansion for this node.
      node->expansion.finalize();

      hpx_lco_set(parms->partdone, 0, nullptr, HPX_NULL, HPX_NULL);
      return HPX_SUCCESS;
    }

    // NOTE: We make use of the SMP only form here. This will not work in
    // distributed.

    // partition sources
    source_t *source_parts{nullptr};
    assert(hpx_gas_try_pin(parms->sources.data(), (void **)&source_parts));

    Source *splits[9] { };
    splits[0] = source_parts;
    splits[8] = &source_parts[parms->sources.n()];

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

    hpx_gas_unpin(parms->sources.data());

    // Find some counts
    sourceref_t cparts[8] { };
    int n_children{0};
    {
      int n_offset{0};
      for (int i = 0; i < 8; ++i) {
        int n_per_child = splits[i + 1] - splits[i];
        if (n_per_child) {
          ++n_children;
          cparts[i] = parms->sources.slice(n_offset, n_per_child);
        }
        n_offset += n_per_child;
      }
    }

    hpx_addr_t childpartdone = hpx_lco_and_new(n_children);
    assert(childpartdone != HPX_NULL);

    // Many arguments to the children's partition are the same
    PartitionParams args{ };
    args.partdone = childpartdone;
    args.limit = parms->limit;
    args.n_digits = parms->n_digits;
    for (int i = 0; i < 8; ++i) {
      if (!cparts[i].valid()) {
        node->child[i] = HPX_NULL;
        continue;
      }

      Index cidx{node->idx.child(i)};
      sourcenode_t thisnode{hpx_thread_current_target()};
      sourcenode_t kid{node->root_geo, cidx, node->method, &thisnode};
      node->child[i] = kid.data();

      // We only set the stuff that changes here
      args.sources = cparts[i];
      hpx_call(kid.data(), partition_, HPX_NULL, &args, sizeof(args));
    }

    int expdigits = parms->n_digits;
    hpx_call_when(childpartdone, hpx_thread_current_target(),
                  child_done_, parms->partdone, &childpartdone, &expdigits);

    return HPX_SUCCESS;
  }

  // The object member data
  mutable Data *local_;
  hpx_addr_t data_;

  // The actions for this class
  static hpx_action_t node_delete_;
  static hpx_action_t self_delete_;
  static hpx_action_t child_done_;
  static hpx_action_t partition_;
};


template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t SourceNode<S, T, E, M>::node_delete_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t SourceNode<S, T, E, M>::self_delete_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t SourceNode<S, T, E, M>::child_done_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t SourceNode<S, T, E, M>::partition_ = HPX_ACTION_NULL;


} // namespace dashmm


#endif // __DASHMM_SOURCE_NODE_H__
