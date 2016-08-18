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


#ifndef __DASHMM_EXPANSION_LCO_H__
#define __DASHMM_EXPANSION_LCO_H__


/// \file include/dashmm/expansionlco.h
/// \brief Interface to Expansion LCO


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
#include "dashmm/shareddata.h"
#include "dashmm/targetlco.h"
#include "dashmm/types.h"
#include "dashmm/viewset.h"


namespace dashmm {


/// Forward declaration of Evaluator so that we can become friends
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class Method,
          typename DistroPolicy>
class Evaluator;


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
/// The referred to object is actually a user-defined LCO that manages the
/// potentially concurrent contribution to the expansion.
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class Method,
          typename DistroPolicy>
class ExpansionLCO {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Method<Source, Target, Expansion, DistroPolicy>;

  using sourceref_t = ArrayRef<Source>;
  using targetref_t = ArrayRef<Target>;
  using targetlco_t = TargetLCO<Source, Target, Expansion, Method,
                                DistroPolicy>;

  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, Method,
                                      DistroPolicy>;

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
  /// \param domain - the domain geometry
  /// \param index - the index of the node containing this expansion
  /// \param expand - initial data for the expansion object
  /// \param where - an addess at the locality where the expansion LCO should be
  ///                created
  ExpansionLCO(int n_in, int n_out, SharedData<DomainGeometry> domain,
               Index index, std::unique_ptr<expansion_t> expand,
               hpx_addr_t where) {
    assert(expand != nullptr);

    ViewSet views = expand->get_all_views();
    size_t bytes = views.bytes();
    size_t less_edges = sizeof(Header) + bytes;

    Header *input_data = reinterpret_cast<Header *>(new char[less_edges]);
    input_data->yet_to_arrive = n_in + 1; // to account for setting out edges
    input_data->expansion_size = bytes;
    input_data->domain = domain;
    input_data->index = index;
    input_data->out_edge_count = n_out;

    WriteBuffer inbuf{input_data->payload, bytes};
    views.serialize(inbuf);

    hpx_addr_t retval{HPX_NULL};
    hpx_call_sync(where, create_from_expansion_, &retval, sizeof(retval),
                  input_data, less_edges);
    delete [] input_data;

    data_ = retval;

    // setup the out edge action
    assert(data_ != HPX_NULL);
    if (n_out != 0) {
      int unused = 0; 
      hpx_call_when(data_, data_, spawn_out_edges_, HPX_NULL, &unused);
    }
  }

  /// Destroy the GAS data referred by the object.
  void destroy() {
    if (data_ != HPX_NULL) {
      hpx_lco_delete_sync(data_);
      data_ = HPX_NULL;
    }
  }

  /// Return the global address of the referred data.
  hpx_addr_t lco() const {return data_;}

  /// Is the object currently referring to global data?
  bool valid() const {return data_ != HPX_NULL;}

  /// Accuracy of expansion
  //int accuracy() const {return 0; } //n_digits_;} //REMOVE

  // TODO: make this more direct. There is no reason to send a parcel since
  // we are going to be local.
  /// Set this expansion with the multipole expansion of the given sources
  ///
  /// This will set the expansion with the multipole expansion computed for
  /// the given @p sources, and with the given @p center. Note that this is
  /// an asynchronous operation. The contribution will be scheduled, but will
  /// not necessarily be complete when this function returns.
  ///
  /// The expansion resulting from the provided sources is added to the
  /// expansion represented by this object.
  ///
  /// \param center - the center of the computed expansion
  /// \param sources - the source data
  /// \param n_src - the number of sources
  void S_to_M(Point center, Source *sources, size_t n_src, Index idx) {
    double scale = expansion_t::compute_scale(idx);
    ViewSet views{kNoRoleNeeded, Point{0.0, 0.0, 0.0}, scale}; 
    expansion_t local{views};
    auto multi = local.S_to_M(center, sources, &sources[n_src]);
    contribute(std::move(multi));
  }

  // TODO: make this more direct. There is no reason to send a parcel since
  // we are going to be local.
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
  void S_to_L(Point center, Source *sources, size_t n_src, Index idx) {
    double scale = expansion_t::compute_scale(idx);     
    ViewSet views{kNoRoleNeeded, Point{0.0, 0.0, 0.0}, scale}; 
    expansion_t local{views};
    auto multi = local.S_to_L(center, sources, &sources[n_src]);
    contribute(std::move(multi));
  }

  // TODO: make this more direct. There is no reason to send a parcel since
  // we are going to be local.
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
  /// \param bytes - the size of the input serialized expansion
  /// \param payload - the serialized expansion data
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
  /// \param ops - the operations
  /// \param targets - the target LCO addresses
  void set_out_edge_data(std::vector<DAGEdge> &edges) {
    int n_out = edges.size();
    size_t bytes = sizeof(OutEdgeRecord) * n_out + sizeof(int) * 2;
    char *input_data = new char[bytes];
    assert(input_data);

    int *codes = reinterpret_cast<int *>(input_data);
    codes[0] = SetOpCodes::kOutEdges;
    codes[1] = n_out;

    if (n_out) {
      std::sort(edges.begin(), edges.end(), DAG::compare_edge_locality);

      OutEdgeRecord *records =
          reinterpret_cast<OutEdgeRecord *>(input_data + sizeof(int) * 2);

      for (int i = 0; i < n_out; ++i) {
        records[i].op = edges[i].op;
        records[i].target = edges[i].target->global_addx;
        records[i].tidx = edges[i].target->idx;
        records[i].locality = edges[i].target->locality;
      }
    }

    hpx_lco_set_lsync(data_, bytes, input_data, HPX_NULL);

    delete [] input_data;
  }

  /// Result the underlying LCO
  ///
  /// This will not only reset the underlying LCO, but will also perform an
  /// 'empty' set of the out edge data.
  void reset() {
    if (data_ == HPX_NULL) return;

    hpx_lco_reset_sync(data_);

    std::vector<DAGEdge> empty{};
    set_out_edge_data(empty);
  }

 private:
  // Give the evaluator access so that it might register our actions
  friend class Evaluator<Source, Target, Expansion, Method, DistroPolicy>;

  ///////////////////////////////////////////////////////////////////
  // Types used internally
  ///////////////////////////////////////////////////////////////////

  /// Part of the internal representation of the Expansion LCO
  ///
  /// Expansion are user-defined LCOs. The data they contain are this object
  /// and the serialized expansion.
  struct Header {
    int yet_to_arrive;
    size_t expansion_size;
    SharedData<DomainGeometry> domain;
    Index index;
    int out_edge_count;
    char payload[];
  };

  /// This stores the data needed to serve the out edge once the LCO is set
  struct OutEdgeRecord {
    Operation op;
    hpx_addr_t target;
    Index tidx;
    int locality;
  };

  /// Operation codes for the LCOs set operation
  enum SetOpCodes {
    kContribute,
    kOutEdges
  };


  ///////////////////////////////////////////////////////////////////
  // LCO Implementation
  ///////////////////////////////////////////////////////////////////

  /// Initialization handler for Expansion LCOs
  ///
  /// The input header is copied directly into the LCO's header. This copies
  /// the metadata as well as the payload (the initial value of the expansion
  /// and the out edges).
  static void init_handler(Header *head, size_t bytes,
                           Header *init, size_t init_bytes) {
    assert(bytes == init_bytes + sizeof(OutEdgeRecord) * init->out_edge_count);
    memcpy(head, init, init_bytes);
  }

  /// The set operation handler for the Expansion LCO
  ///
  /// Set will either save the out edges, or add the input expansion to
  /// the expansion stored in this LCO.
  static void operation_handler(Header *lhs, void *rhs, size_t bytes) {
    ReadBuffer input{static_cast<char *>(rhs), bytes};

    int *code = input.interpret<int>();

    // decrement the counter
    lhs->yet_to_arrive -= 1;
    assert(lhs->yet_to_arrive >= 0);

    switch (*code) {
      case SetOpCodes::kContribute:
        {
          ViewSet views{};
          views.interpret(input);
          expansion_t incoming{views};

          // create the expansion from the payload
          ReadBuffer here{lhs->payload, lhs->expansion_size};
          ViewSet here_views{};
          here_views.interpret(here);
          expansion_t expand{here_views};

          // add the one to the other
          expand.add_expansion(&incoming);

          // release the data, because these objects do not actually own it
          expand.release();
          incoming.release();
        }
        break;
      case SetOpCodes::kOutEdges:
        {
          int n_edges{};
          input.read(&n_edges);

          if (n_edges) {
            // Usage Check
            assert(sizeof(OutEdgeRecord) * n_edges + sizeof(int) * 2 == bytes);
            // We use out edge count here for the size so that an incorrect use
            // does not write outside of bounds.
            WriteBuffer dest{&lhs->payload[lhs->expansion_size],
                             sizeof(OutEdgeRecord) * lhs->out_edge_count};
            dest.write(input);
          }
        }

        break;
    }
  }

  /// The predicate to detect triggering of the Expansion LCO
  ///
  /// If all contributions have arrived, we can now trigger.
  static bool predicate_handler(Header *i, size_t bytes) {
    return (i->yet_to_arrive == 0);
  }


  ///////////////////////////////////////////////////////////////////
  // Other related actions
  ///////////////////////////////////////////////////////////////////

  static int spawn_out_edges_handler(int n_digits) {
    hpx_addr_t lco_ = hpx_thread_current_target();
    // HACK: This action is local to the expansion, so we getref here with
    // whatever as the size and things are okay...
    Header *head{nullptr};
    hpx_lco_getref(lco_, 1, (void **)&head);

    // Shortcut to the work in the case of a single locality
    if (hpx_get_num_ranks() == 1) {
      spawn_out_edges_work(head);
      hpx_lco_release(lco_, head);
      return HPX_SUCCESS;
    }

    // Make a scratch space for the sends
    size_t total_size = sizeof(Header) + head->expansion_size +
                        sizeof(OutEdgeRecord) * head->out_edge_count;
    // TODO increase by sizeof(hpx_addr_t) when return messages happen
    // Actually, this will need to be two of these. One for the tree
    // local pointer array, and one for the return address - actually, if
    // we do the return as a continuation action, we do not need two
    Header *scratch = reinterpret_cast<Header *>(new char [total_size]);
    memcpy(scratch, head, total_size);

    // Once we have the copy, we are done with the LCO data
    hpx_lco_release(lco_, head);

    // Get reference to out edges
    OutEdgeRecord *out_edges = reinterpret_cast<OutEdgeRecord *>(
        &scratch->payload[scratch->expansion_size]);

    // make a local copy of the edge count
    int save_count = scratch->out_edge_count;

    // loop over the sorted edges
    int my_rank = hpx_get_my_rank();
    OutEdgeRecord *begin = out_edges;
    OutEdgeRecord *end = &out_edges[save_count];
    while (begin != end) {
      int curr_rank = begin->locality;

      // Find end of current locality's edges
      OutEdgeRecord *curr = begin;
      while (curr != end && curr->locality == curr_rank) {
        curr += 1;    // look at next record
      }

      //move into start of range
      scratch->out_edge_count = curr - begin;
      size_t send_edges_size = sizeof(OutEdgeRecord) * scratch->out_edge_count;
      if (begin != out_edges) {
        std::move(begin, curr, out_edges);
      }

      // Now send the parcel or do the work
      if (curr_rank == my_rank) {
        spawn_out_edges_work(scratch); 
      } else {
        size_t message_size = sizeof(Header) + scratch->expansion_size
                                           + send_edges_size;
                                    // TODO + (1/2) x hpx_addr_t eventually

        hpx_parcel_t *parc = hpx_parcel_acquire(scratch, message_size);
        hpx_parcel_set_action(parc, spawn_out_edges_from_remote_);
        hpx_parcel_set_target(parc, HPX_THERE(curr_rank));
        // TODO eventually we may need a continuation action and target

        // NOTE: we need to wait for local completion because we are going to
        // modify the buffer in place for the next locality
        hpx_parcel_send_sync(parc);
      }

      //advance
      begin = curr;
    }

    delete [] scratch;

    // done
    return HPX_SUCCESS;
  }

  static int spawn_out_edges_from_remote_handler(Header *head,
                                                 size_t msg_size) {
    // Detect if the edges have unknown target addresses and lookup the
    // correct edges
    // TODO: this is currently not something that will need to happen

    // Send the response message if any lookups were performed
    // TODO: this is not needed just yet. NOTE: the return message will need
    // to be sure that the other sends have finished... Is this an argument
    // for copies on the origin side?
    // TODO: decide if this is a continuation action, or a return action
    // sent from here

    spawn_out_edges_work(head); 

    return HPX_SUCCESS;
  }

  static void spawn_out_edges_work(Header *head) {
    ViewSet views{};
    if (head->out_edge_count > 0) {
      ReadBuffer inbuf{(char *)head->payload, head->expansion_size};
      views.interpret(inbuf);
    }

    // Loop over the out edges, and spawn the work
    OutEdgeRecord *out_edges =
      reinterpret_cast<OutEdgeRecord *>(&head->payload[head->expansion_size]);

    // TODO: should we sort so that the ->T operations go last?
    // Or do two loops, one that skips the ->T and one that does only ->T
    // Explore this when we can measure performance more accurately, or with
    // more detail.
    for (int i = 0; i < head->out_edge_count; ++i) {
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

  static void m_to_m_out_edge(Header *head, const ViewSet &views,
                              hpx_addr_t target) {
    LocalData<DomainGeometry> geo = head->domain.value();
    double s_size = geo->size_from_level(head->index.level());
    int from_child = head->index.which_child();

    expansion_t lexp{views};
    auto translated = lexp.M_to_M(from_child, s_size);
    lexp.release();

    expansionlco_t destination{target};
    destination.contribute(std::move(translated));
  }

  static void m_to_l_out_edge(Header *head, const ViewSet &views,
                              hpx_addr_t target, Index tidx) {
    LocalData<DomainGeometry> geo = head->domain.value();
    double s_size = geo->size_from_level(head->index.level());

    // translate the source expansion
    expansion_t lexp{views};
    auto translated = lexp.M_to_L(head->index, s_size, tidx);
    lexp.release();

    expansionlco_t lco{target};
    lco.contribute(std::move(translated));
  }

  static void l_to_l_out_edge(Header *head, const ViewSet &views,
                              hpx_addr_t target, Index tidx) {
    int to_child = tidx.which_child();
    LocalData<DomainGeometry> geo = head->domain.value();
    double t_size = geo->size_from_level(tidx.level());

    expansion_t lexp{views};
    auto translated = lexp.L_to_L(to_child, t_size);
    lexp.release();

    expansionlco_t total{target};
    total.contribute(std::move(translated));
  }

  static void m_to_t_out_edge(Header *head, hpx_addr_t target) {
    // NOTE: we do not put in the correct number of targets. This is fine
    // because contribute_M_to_T does not rely on this information.
    targetlco_t destination{target, 0}; 
    destination.contribute_M_to_T(head->expansion_size, head->payload);  
  }

  static void l_to_t_out_edge(Header *head, hpx_addr_t target) {
    // NOTE: we do not put in the correct number of targets. This is fine
    // because contribute_L_to_T does not rely on this information.
    targetlco_t destination{target, 0};
    destination.contribute_L_to_T(head->expansion_size, head->payload); 
  }

  static void m_to_i_out_edge(Header *head, const ViewSet &views,
                              hpx_addr_t target) {
    expansion_t lexp{views};
    auto translated = lexp.M_to_I(head->index);
    lexp.release();

    expansionlco_t lco{target};
    lco.contribute(std::move(translated));
  }

  static void i_to_i_out_edge(Header *head, const ViewSet &views,
                              hpx_addr_t target, Index tidx) {
    LocalData<DomainGeometry> geo = head->domain.value();
    double s_size = geo->size_from_level(head->index.level());

    expansion_t lexp{views};
    auto translated = lexp.I_to_I(head->index, s_size, tidx);
    lexp.release();

    expansionlco_t lco{target};
    lco.contribute(std::move(translated));
  }

  static void i_to_l_out_edge(Header *head, const ViewSet &views,
                              hpx_addr_t target, Index tidx) {
    LocalData<DomainGeometry> geo = head->domain.value();
    double t_size = geo->size_from_level(tidx.level());

    expansion_t lexp{views};
    auto translated = lexp.I_to_L(tidx, t_size);
    lexp.release();

    expansionlco_t lco{target};
    lco.contribute(std::move(translated));
  }

  static int create_from_expansion_handler(Header *payload, size_t bytes) {
    size_t total = bytes + sizeof(OutEdgeRecord) * payload->out_edge_count;
    hpx_addr_t gdata = hpx_lco_user_new(total, init_, operation_,
                                        predicate_, payload, bytes);
    assert(gdata != HPX_NULL);
    return HPX_THREAD_CONTINUE(gdata);
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
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t ExpansionLCO<S, T, E, M, D>::init_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t ExpansionLCO<S, T, E, M, D>::operation_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t ExpansionLCO<S, T, E, M, D>::predicate_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t ExpansionLCO<S, T, E, M, D>::spawn_out_edges_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t ExpansionLCO<S, T, E, M, D>::spawn_out_edges_from_remote_ =
    HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t ExpansionLCO<S, T, E, M, D>::create_from_expansion_ =
    HPX_ACTION_NULL;


} // namespace dashmm


#endif // __DASHMM_EXPANSION_REF_H__
