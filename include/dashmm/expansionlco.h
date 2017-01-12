// =============================================================================
//  This file is part of:
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  DASHMM is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  DASHMM is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with DASHMM. If not, see <http://www.gnu.org/licenses/>.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================

#ifndef __DASHMM_EXPANSION_LCO_H__
#define __DASHMM_EXPANSION_LCO_H__


/// \file
/// \brief Interface to Expansion LCO


#include <cstring>

#include <algorithm>
#include <memory>
#include <vector>

#include <hpx/hpx.h>

#include "dashmm/arrayref.h"
#include "dashmm/buffer.h"
#include "dashmm/dag.h"
#include "dashmm/domaingeometry.h"
#include "dashmm/index.h"
#include "dashmm/point.h"
#include "dashmm/rankwise.h"
#include "dashmm/shareddata.h"
#include "dashmm/targetlco.h"
#include "dashmm/traceevents.h"
#include "dashmm/types.h"
#include "dashmm/viewset.h"


namespace dashmm {


template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class> class Method>
class DualTree;


/// Forward declaration of registrar so that we can become friends
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class> class Method>
class ExpansionLCORegistrar;


/// Expansion LCO
///
/// This object is a thin wrapper around the address of a user-defined LCO
/// whose data is the serialized form of an expansion. The main purpose of
/// this object is to hide HPX-5 from end-users that are not interested in the
/// details of HPX-5.
///
/// This object attempts to have as many of the same methods as the
/// underlying Expansion object. However, some of the function signatures are
/// different, and some are missing.
///
/// The referred to object is, in the HPX-5 parlance, a user-defined LCO that
/// manages the potentially concurrent contribution to the expansion.
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class> class Method>
class ExpansionLCO {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Method<Source, Target, Expansion>;

  using targetlco_t = TargetLCO<Source, Target, Expansion, Method>;
  using dualtree_t = DualTree<Source, Target, Expansion, Method>;

  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, Method>;


  /// Construct the expansion from a given global address.
  ExpansionLCO(hpx_addr_t addr) : data_{addr} {}

  /// Construct the expansion LCO from expansion data
  ///
  /// This will create the expansion LCO at the locality specified by where,
  /// with the provided data. The number of inputs will be increased by one
  /// to account for the fact that the out edges will need to be set after the
  /// expansion LCOs are created.
  ///
  /// \param n_in - the number of input edges to the expansion
  /// \param n_out - the number of output edges from the expansion
  /// \param index - the index of the node containing this expansion
  /// \param expand - initial data for the expansion object
  /// \param where - an addess at the locality where the expansion LCO should be
  ///                created
  /// \param rwtree - the global address of the dual tree
  ExpansionLCO(int n_in, int n_out,
               Index index, std::unique_ptr<expansion_t> expand,
               hpx_addr_t rwtree) {
    assert(expand != nullptr);

    ViewSet views = expand->get_all_views();
    size_t bytes = views.bytes();

    data_ = hpx_lco_user_new(sizeof(Header), init_, operation_, predicate_,
                             &bytes, sizeof(bytes));
    assert(data_ != HPX_NULL);

    void *lva{nullptr};
    assert(hpx_gas_try_pin(data_, &lva));
    Header *ldata = static_cast<Header *>(hpx_lco_user_get_user_data(lva));

    ldata->yet_to_arrive = n_in + 1; // to account for setting out edges
    ldata->index = index;
    ldata->rwaddr = rwtree;
    ldata->out_edge_count = n_out;
    ldata->data = new char[bytes + sizeof(OutEdgeRecord) * n_out];

    WriteBuffer inbuf{ldata->data, bytes};
    views.serialize(inbuf);

    hpx_gas_unpin(data_);

    // setup the out edge action
    if (n_out != 0) {
      int unused = 0;
      hpx_call_when(data_, data_, spawn_out_edges_, HPX_NULL, &unused);
    }
  }

  /// Destroy the GAS data referred by the object.
  void destroy() {
    if (data_ != HPX_NULL) {
      void *lva{nullptr};
      assert(hpx_gas_try_pin(data_, &lva));
      Header *ldata = static_cast<Header *>(hpx_lco_user_get_user_data(lva));
      delete [] ldata->data;
      hpx_gas_unpin(data_);

      hpx_lco_delete_sync(data_);
      data_ = HPX_NULL;
    }
  }

  /// Return the global address of the referred data.
  hpx_addr_t lco() const {return data_;}

  /// Is the object currently referring to global data?
  bool valid() const {return data_ != HPX_NULL;}

  /// Set this expansion with the multipole expansion of the given sources
  ///
  /// This will set the expansion with the multipole expansion computed for the
  /// given @p sources. As this M only expects 1 input, the implementation here
  /// directly operates on the memory reserved for the expansion data inside
  /// this object. The set operation performed simply decrements the internal
  /// counter of the LCO. 
  ///
  /// \param sources - the source data
  /// \param n_src - the number of sources
  void S_to_M(Source *sources, size_t n_src) {
    EVENT_TRACE_DASHMM_STOM_BEGIN();

    void *lva{nullptr}; 
    assert(hpx_gas_try_pin(data_, &lva)); 
    Header *head = static_cast<Header *>(hpx_lco_user_get_user_data(lva)); 
    ReadBuffer here{head->data, head->expansion_size}; 
    ViewSet views{}; 
    views.interpret(here); 
    expansion_t local{views}; 
    local.S_to_M(sources, &sources[n_src]); 
    local.release(); 
    hpx_gas_unpin(data_); 
    int code = SetOpCodes::kStoM; 
    hpx_lco_set_lsync(data_, sizeof(int), &code, HPX_NULL); 

    EVENT_TRACE_DASHMM_STOM_END();
  }

  /// Set this expansion with the local expansion of the given sources
  ///
  /// This will set the expansion with the local expansion computed for
  /// the given @p sources, and with the given @p center. Note that this is
  /// an asynchronous operation. The contribution will be scheduled, but will
  /// not necessarily be complete when this function returns.
  ///
  /// Note that this does not set the expansion to be equal to the computed
  /// local expansion. Instead, it adds the computed local to the
  /// current contents of the expansion.
  ///
  /// \param center - the center of the computed expansion
  /// \param sources - the source data
  /// \param n_src - the number of sources
  /// \param idx - the index of the node for which this is an expansion
  void S_to_L(Point center, Source *sources, size_t n_src, Index idx) {
    EVENT_TRACE_DASHMM_STOL_BEGIN();
    double scale = expansion_t::compute_scale(idx);
    ViewSet views{kNoRoleNeeded, Point{0.0, 0.0, 0.0}, scale};
    expansion_t local{views};
    auto multi = local.S_to_L(center, sources, &sources[n_src]);
    contribute(std::move(multi));
    EVENT_TRACE_DASHMM_STOL_END();
  }

  /// Apply effect of sources to targets
  ///
  /// This will compute the effect of the given @p sources on the given
  /// @p targets. Note that this is an asynchronous operation. This may
  /// return before the contribution to the targets has been computed.
  ///
  /// \param sources - the source data
  /// \param n_sources - the number of sources
  /// \param targets - a reference to the target points
  void S_to_T(Source *sources, size_t n_sources, targetlco_t targets) const {
    targets.contribute_S_to_T(n_sources, sources);
  }

  /// Contribute to the referred expansion
  ///
  /// This will call the appropriate set operation on the referred LCO. This
  /// will result in the add_expansion method of the expansion being called.
  ///
  /// \param expand - the expansion to contribute
  void contribute(std::unique_ptr<expansion_t> &&expand) {
    ViewSet views = expand->get_all_views();
    size_t bytes = views.bytes();

    hpx_parcel_t *parc = hpx_parcel_acquire(nullptr, bytes + sizeof(int));
    assert(parc != nullptr);

    hpx_parcel_set_action(parc, hpx_lco_set_action);
    hpx_parcel_set_target(parc, data_);

    WriteBuffer parcbuf{(char *)hpx_parcel_get_data(parc), bytes + sizeof(int)};

    int opcode = SetOpCodes::kContribute;
    parcbuf.write((char *)&opcode, sizeof(opcode));
    views.serialize(parcbuf);

    // We do not need local completion because we do not own this parcel, so
    // we will not delete it. And so since we are already copied into the
    // parcel, we can go ahead and move on with life.
    hpx_parcel_send(parc, HPX_NULL);

    // NOTE: No release is needed. The argument will go out of scope at this
    // point, and the data will be freed.
  }

  /// Set the out edge data for this expansion LCO
  ///
  /// This will set the underlying LCO with the correct out edge data.
  ///
  /// \param edges - the out edge data
  void set_out_edge_data(std::vector<DAGEdge> &edges) {
    void *lva{nullptr};
    assert(hpx_gas_try_pin(data_, &lva));
    Header *ldata = static_cast<Header *>(hpx_lco_user_get_user_data(lva));

    int n_out = edges.size();
    assert(ldata->out_edge_count == n_out);

    if (n_out) {
      std::sort(edges.begin(), edges.end(), DAG::compare_edge_locality);

      OutEdgeRecord *records =
        reinterpret_cast<OutEdgeRecord *>(ldata->data + ldata->expansion_size);

      for (int i = 0; i < n_out; ++i) {
        records[i].op = edges[i].op;
        records[i].target = edges[i].target->global_addx;
        records[i].tidx = edges[i].target->idx;
        records[i].locality = edges[i].target->locality;
      }
    }

    hpx_gas_unpin(data_);
    int code = SetOpCodes::kOutEdges;
    hpx_lco_set_lsync(data_, sizeof(int), &code, HPX_NULL);
  }

  /// Reset the underlying LCO
  ///
  /// This will not only reset the underlying LCO, but will also perform an
  /// 'empty' set of the out edge data. The intent of this method is for use
  /// in cases where a DAG is reused multiple times. At the moment, this
  /// functionality is not available in DASHMM, but will be soon.
  void reset() {
    if (data_ == HPX_NULL) return;

    hpx_lco_reset_sync(data_);

    std::vector<DAGEdge> empty{};
    set_out_edge_data(empty);
  }

 private:
  // Give the registrar access so that it might register our actions
  friend class ExpansionLCORegistrar<Source, Target, Expansion, Method>;

  ///////////////////////////////////////////////////////////////////
  // Types used internally
  ///////////////////////////////////////////////////////////////////

  /// This stores the data needed to serve the out edge once the LCO is set
  struct OutEdgeRecord {
    Operation op;
    hpx_addr_t target;
    Index tidx;
    int locality;
  };

  /// Part of the internal representation of the Expansion LCO
  ///
  /// Expansion are user-defined LCOs. The data they contain are this object
  /// and the serialized expansion. The expansion_data points to the serialized
  /// form of the expansion. This is done through the ViewSet object, and
  /// details on the exact format can be found with the ViewSet documentation.
  struct Header {
    int yet_to_arrive;
    int out_edge_count;
    size_t expansion_size;
    Index index;
    hpx_addr_t rwaddr;
    char *data;
  };

  /// Operation codes for the LCOs set operation
  enum SetOpCodes {
    kContribute,
    kStoM, 
    kOutEdges
  };


  ///////////////////////////////////////////////////////////////////
  // LCO Implementation
  ///////////////////////////////////////////////////////////////////

  /// Initialization handler for Expansion LCOs
  ///
  /// This updates the expansion size, and also zeros out the expansion
  /// data if the pointer is valid.
  ///
  /// \param head - the address of the LCO data
  /// \param bytes - the size of the LCO data
  /// \param init - the initialization data for the LCO
  /// \param init_bytes - the size of the initialization data for the LCO
  static void init_handler(Header *head, size_t bytes,
                           size_t *init, size_t init_bytes) {
    head->expansion_size = *init;
  }

  /// The set operation handler for the Expansion LCO
  ///
  /// Set will either add the input expansion to the expansion referenced in
  /// this LCO, or simply decrement the counter that monitors the status of the
  /// LCO
  ///
  /// In both cases, @p rhs begins with an integer indicating which sort of set
  /// was called. If the set is to accumulate expansion, then there is a
  /// serialized ViewSet following the integer code. This buffer is deserialized
  /// into an expansion, and then added to the expansion referenced by this
  /// LCO.
  ///
  /// \param lhs - the address of this LCO's data
  /// \param rhs - the input buffer
  /// \param bytes - the size of the input buffer
  static void operation_handler(Header *lhs, void *rhs, size_t bytes) {
    ReadBuffer input{static_cast<char *>(rhs), bytes};

    int *code = input.interpret<int>();

    // decrement the counter
    lhs->yet_to_arrive -= 1;
    assert(lhs->yet_to_arrive >= 0);

    switch (*code) {
    case SetOpCodes::kContribute:
      {
        EVENT_TRACE_DASHMM_ELCO_BEGIN();
        ViewSet views{};
        views.interpret(input);
        expansion_t incoming{views};
        
        ReadBuffer here{lhs->data, lhs->expansion_size};
        ViewSet here_views{};
        here_views.interpret(here);
        expansion_t expand{here_views};
        
        // add the one to the other
        expand.add_expansion(&incoming);
        
        // release the data, because these objects do not actually own it
        expand.release();
        incoming.release();
        EVENT_TRACE_DASHMM_ELCO_END();
      }
      break;
    case SetOpCodes::kStoM: 
      break;
    case SetOpCodes::kOutEdges:
      break;
    }
  }

  /// The predicate to detect triggering of the Expansion LCO
  ///
  /// If all contributions have arrived, and the out edge data has been set,
  /// we can trigger.
  static bool predicate_handler(Header *i, size_t bytes) {
    return (i->yet_to_arrive == 0);
  }


  ///////////////////////////////////////////////////////////////////
  // Other related actions
  ///////////////////////////////////////////////////////////////////

  /// Spawn the work at the out edges of this LCO
  ///
  /// Once the LCO is triggered, it will perform the actions required by the
  /// out edges of that LCO. This happens typically in two steps. The first is
  /// bundling of the out edges that go to LCOs at the same locality together
  /// and sending the data across the network a single time. Then at the
  /// remote side, the addresses are looked up and the work is performed.
  /// See spawn_out_edges_from_remote_handler and spawn_out_edges_work for
  /// more details.
  ///
  /// \param unused - is not used
  ///
  /// \returns - HPX_SUCCESS
  static int spawn_out_edges_handler(int unused) {
    hpx_addr_t lco_ = hpx_thread_current_target();
    // HACK: This action is local to the expansion, so we getref here with
    // whatever as the size and things are okay...
    Header *head{nullptr};
    hpx_lco_getref(lco_, 1, (void **)&head);

    // Shortcut to the work in the case of a single locality
    if (hpx_get_num_ranks() == 1) {
      spawn_out_edges_work(head, 0, head->out_edge_count - 1);
      hpx_lco_release(lco_, head);
      return HPX_SUCCESS;
    }

    // Make a scratch space for the sends
    size_t edgeless = sizeof(Header) + head->expansion_size;
    int out_edge_count = head->out_edge_count;
    size_t edge_size = sizeof(OutEdgeRecord) * out_edge_count;
    size_t total = edgeless + edge_size;

    char *temp = new char[total];
    Header *scratch = reinterpret_cast<Header *>(temp);

    // Fill in the header
    memcpy(scratch, head, sizeof(Header));

    // Fill in the expansion data
    memcpy(temp + sizeof(Header), head->data, head->expansion_size);

    // Address for storing out edges
    OutEdgeRecord *scratch_edges =
        reinterpret_cast<OutEdgeRecord *>(temp + edgeless);

    OutEdgeRecord *out_edges =
        reinterpret_cast<OutEdgeRecord *>(head->data + head->expansion_size);

    // Loop over the sorted edges
    int my_rank = hpx_get_my_rank();
    int begin = 0;

    while (begin != out_edge_count) {
      int curr_rank = out_edges[begin].locality;

      int curr = begin;
      while (curr != out_edge_count &&
             out_edges[curr].locality == curr_rank) {
        curr++; // look at next record
      }

      if (curr_rank == my_rank) {
        // Short cut to do work
        spawn_out_edges_work(head, begin, curr - 1);
      } else {
        int curr_rank_out_edge_count = curr - begin;

        // Update the edge count information
        scratch->out_edge_count = curr_rank_out_edge_count;

        // Copy the edge information into scratch
        memcpy(scratch_edges, &out_edges[begin],
               sizeof(OutEdgeRecord) * curr_rank_out_edge_count);

        // Prepare parcel
        size_t message_size = edgeless +
          sizeof(OutEdgeRecord) * curr_rank_out_edge_count;
        // TODO + (1/2) x hpx_addr_t eventually

        hpx_parcel_t *parc = hpx_parcel_acquire(temp, message_size);
        hpx_parcel_set_action(parc, spawn_out_edges_from_remote_);
        hpx_parcel_set_target(parc, HPX_THERE(curr_rank));

        // NOTE: we need to wait for local completion because we are going to
        // modify the buffer in place for the next locality
        hpx_parcel_send_sync(parc);
      }

      // Advance
      begin = curr;
    }

    delete [] temp;

    // done
    return HPX_SUCCESS;
  }

  /// Action to handle incoming edges from a remote
  ///
  /// This action is spawned when out edges from a remote arrive. It will first
  /// look up any LCO address based on the index in the local tree and then it
  /// will perform the work. The possibility that some might not need to be
  /// looked up relies on the future possible reuse of the DAG. This is a
  /// feature that does not currently exist, but which might in the future.
  ///
  /// The message contains an expansion's data, but with only a subset of the
  /// out edge records.
  ///
  /// \param head - the message data from the remote
  /// \param msg_size - the size of the incoming data
  ///
  /// \returns - HPX_SUCCESS
  static int spawn_out_edges_from_remote_handler(Header *head,
                                                 size_t msg_size) {
    // Detect if the edges have unknown target addresses and lookup the
    // correct edges
    RankWise<dualtree_t> global_tree{head->rwaddr};
    auto tree = global_tree.here();

    // \p head->expansion_data and \p head->out_edge_records currently store the
    // addresses of the sender's side, which are invalid on the receiver's
    // side. Compute correct offsets to set these pointers correctly.
    char *temp = reinterpret_cast<char *>(head);
    head->data = temp + sizeof(Header);
    OutEdgeRecord *out_edges =
      reinterpret_cast<OutEdgeRecord *>(temp + sizeof(Header) +
                                        head->expansion_size);

    for (int i = 0; i < head->out_edge_count; ++i) {
      if (out_edges[i].target == HPX_NULL) {
        out_edges[i].target = tree->lookup_lco_addx(out_edges[i].tidx,
                                                    out_edges[i].op);
      }
    }

    spawn_out_edges_work(head, 0, head->out_edge_count - 1);

    return HPX_SUCCESS;
  }

  /// Routine driving the actual work of the incoming edges
  ///
  /// This is called either after the address lookup, or directly for those
  /// edges local to the initiating LCO.
  ///
  /// \param head - the message data
  /// \param first - the first edge to be processed
  /// \param last - the last edge to be proceed
  static void spawn_out_edges_work(Header *head, int first, int last) {
    ViewSet views{};

    if (head->out_edge_count > 0) {
      ReadBuffer inbuf{head->data, head->expansion_size};
      views.interpret(inbuf);
    }

    OutEdgeRecord *out_edges =
        reinterpret_cast<OutEdgeRecord *>(head->data + head->expansion_size);

    for (int i = first; i <= last; ++i) {
      switch(out_edges[i].op) {
        case Operation::MtoM:
          m_to_m_out_edge(head, views, out_edges[i].target);
          break;
        case Operation::MtoL:
          m_to_l_out_edge(head, views, out_edges[i].target, out_edges[i].tidx);
          break;
        case Operation::LtoL:
          l_to_l_out_edge(head, views, out_edges[i].target, out_edges[i].tidx);
          break;
        case Operation::MtoT:
          m_to_t_out_edge(head, out_edges[i].target);
          break;
        case Operation::LtoT:
          l_to_t_out_edge(head, out_edges[i].target);
          break;
        case Operation::MtoI:
          m_to_i_out_edge(head, views, out_edges[i].target);
          break;
        case Operation::ItoI:
          i_to_i_out_edge(head, views, out_edges[i].target, out_edges[i].tidx);
          break;
        case Operation::ItoL:
          i_to_l_out_edge(head, views, out_edges[i].target, out_edges[i].tidx);
          break;
        default:
          assert(0 && "Impossible operation during out edge spawn");
          break;
      }
    }
  }

  /// Serve an M->M edge
  ///
  /// \param head - the incoming data
  /// \param views - ViewSet serving the expansion
  /// \param target - global address of target LCO
  static void m_to_m_out_edge(Header *head, const ViewSet &views,
                              hpx_addr_t target) {
    EVENT_TRACE_DASHMM_MTOM_BEGIN();
    int from_child = head->index.which_child(); 
    expansion_t lexp{views}; 
    auto translated = lexp.M_to_M(from_child); 
    lexp.release();  // Question: why release? (clears the views)

    expansionlco_t destination{target}; 
    destination.contribute(std::move(translated)); 
    EVENT_TRACE_DASHMM_MTOM_END();
  }

  /// Serve an M->L edge
  ///
  /// \param head - the incoming data
  /// \param views - ViewSet serving the expansion
  /// \param target - global address of target LCO
  /// \param tidx - index of target LCO
  static void m_to_l_out_edge(Header *head, const ViewSet &views,
                              hpx_addr_t target, Index tidx) {
    EVENT_TRACE_DASHMM_MTOL_BEGIN();

    // translate the source expansion
    expansion_t lexp{views};
    auto translated = lexp.M_to_L(head->index, tidx);
    lexp.release();

    expansionlco_t lco{target};
    lco.contribute(std::move(translated));
    EVENT_TRACE_DASHMM_MTOL_END();
  }

  /// Serve an L->L edge
  ///
  /// \param head - the incoming data
  /// \param views - ViewSet serving the expansion
  /// \param target - global address of target LCO
  /// \param tidx - index of target LCO
  static void l_to_l_out_edge(Header *head, const ViewSet &views,
                              hpx_addr_t target, Index tidx) {
    EVENT_TRACE_DASHMM_LTOL_BEGIN();
    int to_child = tidx.which_child();

    expansion_t lexp{views};
    auto translated = lexp.L_to_L(to_child); 
    lexp.release();

    expansionlco_t total{target};
    total.contribute(std::move(translated));
    EVENT_TRACE_DASHMM_LTOL_END();
  }

  /// Serve an M->T edge
  ///
  /// \param head - the incoming data
  /// \param target - global address of target LCO
  static void m_to_t_out_edge(Header *head, hpx_addr_t target) {
    // NOTE: we do not put in the correct number of targets. This is fine
    // because contribute_M_to_T does not rely on this information.
    targetlco_t destination{target, 0};
    destination.contribute_M_to_T(head->expansion_size, head->data);
  }

  /// Serve an L->T edge
  ///
  /// \param head - the incoming data
  /// \param target - global address of target LCO
  static void l_to_t_out_edge(Header *head, hpx_addr_t target) {
    // NOTE: we do not put in the correct number of targets. This is fine
    // because contribute_L_to_T does not rely on this information.
    targetlco_t destination{target, 0};
    destination.contribute_L_to_T(head->expansion_size, head->data);
  }

  /// Serve an M->I edge
  ///
  /// \param head - the incoming data
  /// \param views - ViewSet serving the expansion
  /// \param target - global address of target LCO
  static void m_to_i_out_edge(Header *head, const ViewSet &views,
                              hpx_addr_t target) {
    EVENT_TRACE_DASHMM_MTOI_BEGIN();
    expansion_t lexp{views};
    auto translated = lexp.M_to_I(); 
    lexp.release();

    expansionlco_t lco{target};
    lco.contribute(std::move(translated));
    EVENT_TRACE_DASHMM_MTOI_END();
  }

  /// Serve an I->I edge
  ///
  /// \param head - the incoming data
  /// \param views - ViewSet serving the expansion
  /// \param target - global address of target LCO
  /// \param tidx - index of target LCO
  static void i_to_i_out_edge(Header *head, const ViewSet &views,
                              hpx_addr_t target, Index tidx) {
    EVENT_TRACE_DASHMM_ITOI_BEGIN();

    expansion_t lexp{views};
    auto translated = lexp.I_to_I(head->index, tidx);
    lexp.release();

    expansionlco_t lco{target};
    lco.contribute(std::move(translated));
    EVENT_TRACE_DASHMM_ITOI_END();
  }

  /// Serve an I->L edge
  ///
  /// \param head - the incoming data
  /// \param views - ViewSet serving the expansion
  /// \param target - global address of target LCO
  /// \param tidx - index of target LCO
  static void i_to_l_out_edge(Header *head, const ViewSet &views,
                              hpx_addr_t target, Index tidx) {
    EVENT_TRACE_DASHMM_ITOL_BEGIN();

    expansion_t lexp{views};
    auto translated = lexp.I_to_L(tidx); 
    lexp.release();

    expansionlco_t lco{target};
    lco.contribute(std::move(translated));
    EVENT_TRACE_DASHMM_ITOL_END();
  }


  // The functions implementing the user LCO
  static hpx_action_t init_;
  static hpx_action_t operation_;
  static hpx_action_t predicate_;

  // The actions for the various operations
  static hpx_action_t spawn_out_edges_;
  static hpx_action_t spawn_out_edges_from_remote_;
  static hpx_action_t create_from_expansion_;

  hpx_addr_t data_;     // this is the LCO
};

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t ExpansionLCO<S, T, E, M>::init_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t ExpansionLCO<S, T, E, M>::operation_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t ExpansionLCO<S, T, E, M>::predicate_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t ExpansionLCO<S, T, E, M>::spawn_out_edges_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t ExpansionLCO<S, T, E, M>::spawn_out_edges_from_remote_ =
    HPX_ACTION_NULL;


} // namespace dashmm


#endif // __DASHMM_EXPANSION_REF_H__
