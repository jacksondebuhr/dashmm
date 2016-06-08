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


#include <memory>
#include <vector>

#include <hpx/hpx.h>

#include "dashmm/arrayref.h"
#include "dashmm/daginfo.h"
#include "dashmm/domaingeometry.h"
#include "dashmm/index.h"
#include "dashmm/point.h"
#include "dashmm/shareddata.h"
#include "dashmm/targetlco.h"
#include "dashmm/types.h"


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
  ExpansionLCO(hpx_addr_t addr, int n_digits)
      : data_{addr}, n_digits_{n_digits} { }

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

    size_t bytes = expand->bytes();
    size_t total_size = sizeof(Header) + bytes + sizeof(OutEdgeRecord) * n_out;

    Header *input_data = reinterpret_cast<Header *>(new char[total_size]);
    input_data->yet_to_arrive = n_in + 1; // to account for setting out edges
    input_data->expansion_size = bytes;
    input_data->domain = domain;
    input_data->index = index;
    input_data->out_edge_count = n_out;

    int n_digits = expand->accuracy();
    char *ldata = reinterpret_cast<char *>(expand->release());
    memcpy(input_data->payload, ldata, bytes);

    hpx_addr_t retval{HPX_NULL};
    hpx_call_sync(where, create_from_expansion_, &retval, sizeof(retval),
                  input_data, total_size);
    delete [] ldata;

    data_ = retval;
    n_digits_ = n_digits;

    // setup the out edge action
    hpx_call_when(data_, data_, spawn_out_edges_, HPX_NULL, &n_digits);
  }

  /// Destroy the GAS data referred by the object.
  void destroy() {
    if (data_ != HPX_NULL) {
      hpx_lco_delete_sync(data_);
      data_ = HPX_NULL;
      n_digits_ = -1;
    }
  }

  /// Return the global address of the referred data.
  hpx_addr_t data() const {return data_;}

  /// Is the object currently referring to global data?
  bool valid() const {return data_ != HPX_NULL;}

  /// Accuracy of expansion
  int accuracy() const {return n_digits_;}

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
  /// \param sources - a reference to the sources from which to compute the
  ///                  multipole expansion
  /// \param scale - scaling factor
  void S_to_M(Point center, sourceref_t sources, double scale) const {
    int nsrc = sources.n();
    double cx = center.x();
    double cy = center.y();
    double cz = center.z();
    hpx_call(sources.data(), s_to_m_, HPX_NULL, &nsrc, &cx, &cy, &cz,
             &scale, &data_, &n_digits_);
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
  /// \param sources - a reference to the sources from which to compute the
  ///                  local expansion
  /// \param scale - scaling factor
  void S_to_L(Point center, sourceref_t sources, double scale) const {
    int nsrc = sources.n();
    double cx = center.x();
    double cy = center.y();
    double cz = center.z();
    hpx_call(sources.data(), s_to_l_, HPX_NULL, &nsrc, &cx, &cy, &cz, &scale,
             &data_, &n_digits_);
  }

  /// Apply effect of sources to targets
  ///
  /// This will compute the effect of the given @p sources on the given
  /// @p targets. Note that this is an asynchronous operation. This may
  /// return before the contribution to the targets has been computed.
  ///
  /// \param sources - a reference to the source points
  /// \param targets - a reference to the target points
  void S_to_T(sourceref_t sources, targetlco_t targets) const {
    int n_src = sources.n();
    hpx_addr_t tsend = targets.lco();
    int n_trg = targets.n();
    hpx_call(sources.data(), s_to_t_, HPX_NULL, &n_src, &tsend, &n_trg);
  }

  /// Contribute to the referred expansion
  ///
  /// This will call the appropriate set operation on the referred LCO. This
  /// will result in the add_expansion method of the expansion being called.
  ///
  /// \param bytes - the size of the input serialized expansion
  /// \param payload - the serialized expansion data
  void contribute(size_t bytes, char *payload) {
    int *code = reinterpret_cast<int *>(payload);
    code[0] = SetOpCodes::kContribute;
    hpx_lco_set_lsync(data_, bytes, payload, HPX_NULL);
  }

  /// Set the out edge data for this expansion LCO
  ///
  /// This will set the underlying LCO with the correct out edge data.
  ///
  /// \param ops - the operations
  /// \param targets - the target LCO addresses
  void set_out_edge_data(const std::vector<DAGEdge> &edges) {
    int n_out = edges.size();
    size_t bytes = sizeof(OutEdgeRecord) * n_out + sizeof(int) * 2;
    char *input_data = new char[bytes];
    assert(input_data);

    int *codes = reinterpret_cast<int *>(input_data);
    OutEdgeRecord *records =
        reinterpret_cast<OutEdgeRecord *>(input_data + sizeof(int) * 2);

    codes[0] = SetOpCodes::kOutEdges;
    codes[1] = n_out;

    for (int i = 0; i < n_out; ++i) {
      records[i].op = edges[i].op;
      records[i].target = edges[i].target->global_addx;
    }

    hpx_lco_set_lsync(data_, bytes, input_data, HPX_NULL);

    delete [] input_data;
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
  };

  /// Marshalled parameter type for m_to_l_
  struct MtoLParams {
    size_t bytes;
    int n_digits;
    Index index;
    char payload[];
  };

  /// Marshalled parameter type for l_to_l_
  struct LtoLParams {
    size_t bytes;
    int n_digits;
    char payload[];
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
    assert(bytes == init_bytes);
    memcpy(head, init, init_bytes);
  }

  /// The set operation handler for the Expansion LCO
  ///
  /// Set will either save the out edges, or add the input expansion to
  /// the expansion stored in this LCO.
  static void operation_handler(Header *lhs, void *rhs, size_t bytes) {
    int *code = static_cast<int *>(rhs);

    // decrement the counter
    lhs->yet_to_arrive -= 1;
    assert(lhs->yet_to_arrive >= 0);

    switch (code[0]) {
      case SetOpCodes::kContribute:
        {
          expansion_t incoming{rhs, bytes, code[1]};

          // create the expansion from the payload
          int *n_digits = reinterpret_cast<int *>(lhs->payload + sizeof(int));
          expansion_t expand{lhs->payload, lhs->expansion_size, *n_digits};

          // add the one to the other
          expand.add_expansion(&incoming);

          // release the data, because these objects do not actually own it
          expand.release();
          incoming.release();
        }
        break;
      case SetOpCodes::kOutEdges:
        {
          int n_edges = code[1];
          // usage check
          assert(sizeof(OutEdgeRecord) * n_edges + sizeof(int) * 2 == bytes);

          char *dest = &lhs->payload[lhs->expansion_size];
          char *offset = static_cast<char *>(rhs) + sizeof(int) * 2;
          memcpy(dest, offset, sizeof(OutEdgeRecord) * n_edges);
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

  // I think the approach is going to be to do the action immediately,
  // moving through the list in serial. If this ends up being a problem, we
  // shall have to make each into an action, and accept the overhead.
  static int spawn_out_edges_handler(int n_digits) {
    hpx_addr_t lco_ = hpx_thread_current_target();
    // HACK: This action is local to the expansion, so we getref here with
    // whatever as the size and things are okay...
    Header *head{nullptr};
    hpx_lco_getref(lco_, 1, (void **)&head);

    // Loop over the out edges, and spawn the work
    OutEdgeRecord *out_edges =
      reinterpret_cast<OutEdgeRecord *>(&head->payload[head->expansion_size]);
    for (int i = 0; i < head->out_edge_count; ++i) {
      switch(out_edges[i].op) {
        case Operation::MtoM:
          m_to_m_out_edge(head, out_edges[i].target, n_digits);
          break;
        case Operation::MtoL:
          m_to_l_out_edge(head, out_edges[i].target, n_digits);
          break;
        case Operation::LtoL:
          l_to_l_out_edge(head, out_edges[i].target, n_digits);
          break;
        case Operation::MtoT:
          m_to_t_out_edge(head, out_edges[i].target, n_digits);
          break;
        case Operation::LtoT:
          l_to_t_out_edge(head, out_edges[i].target, n_digits);
          break;
        default:
          break;
      }
    }

    hpx_lco_release(lco_, head);

    return HPX_SUCCESS;
  }

  static int s_to_m_handler(Source *sources, int n_src, double cx, double cy,
                            double cz, double scale, hpx_addr_t expand,
                            int n_digits) {
    expansion_t local{nullptr, 0, n_digits};
    std::unique_ptr<expansion_t> multi = std::move(
      local.S_to_M(Point{cx, cy, cx}, sources, &sources[n_src], scale));
    size_t bytes = multi->bytes();
    char *serial = static_cast<char *>(multi->release());

    expansionlco_t total{expand, n_digits};
    total.contribute(bytes, serial);
    delete [] serial;

    return HPX_SUCCESS;
  }

  static int s_to_l_handler(Source *sources, int n_src, double cx, double cy,
                            double cz, double scale, hpx_addr_t expand,
                            int n_digits) {
    expansion_t local{nullptr, 0, n_digits};
    auto multi = local.S_to_L(Point{cx, cy, cz}, sources, &sources[n_src],
                              scale);
    size_t bytes = multi->bytes();
    char *serial = static_cast<char *>(multi->release());

    expansionlco_t total{expand, n_digits};
    total.contribute(bytes, serial);
    delete [] serial;

    return HPX_SUCCESS;
  }

  static void m_to_m_out_edge(Header *head, hpx_addr_t target, int n_digits) {
    LocalData<DomainGeometry> geo = head->domain.value();
    double s_size = geo->size_from_level(head->index.level());
    int from_child = head->index.which_child();

    expansion_t lexp{head->payload, head->expansion_size, n_digits};
    auto translated = lexp.M_to_M(from_child, s_size);
    lexp.release();

    size_t bytes = translated->bytes();
    char *transexpand = reinterpret_cast<char *>(translated->release());

    expansionlco_t destination{target, n_digits};
    destination.contribute(bytes, transexpand);

    delete [] transexpand;
  }

  static void m_to_l_out_edge(Header *head, hpx_addr_t target, int n_digits) {
    // Get a parcel
    size_t msg_size = head->expansion_size + sizeof(MtoLParams);
    hpx_parcel_t *parc = hpx_parcel_acquire(nullptr, msg_size);
    assert(parc != nullptr);

    // Setup the parcel
    hpx_parcel_set_action(parc, m_to_l_);
    hpx_parcel_set_target(parc, target);

    // Fill the parcel payload
    MtoLParams *message = static_cast<MtoLParams *>(hpx_parcel_get_data(parc));
    message->bytes = head->expansion_size;
    message->n_digits = n_digits;
    message->index = head->index;
    memcpy(message->payload, head->payload, head->expansion_size);

    // Send the parcel - this will release the acquired parcel
    hpx_parcel_send_sync(parc);
  }

  static int m_to_l_handler(MtoLParams *parms, size_t UNUSED) {
    hpx_addr_t target = hpx_thread_current_target();

    void *workaround{nullptr};
    assert(hpx_gas_try_pin(target, &workaround));
    Header *trg_data =
        static_cast<Header *>(hpx_lco_user_get_user_data(workaround));
    Index t_index = trg_data->index;
    LocalData<DomainGeometry> geo = trg_data->domain.value();
    hpx_gas_unpin(target);
    double s_size = geo->size_from_level(parms->index.level());

    // translate the source expansion
    expansion_t lexp{parms->payload, parms->bytes, parms->n_digits};
    auto translated = lexp.M_to_L(parms->index, s_size, t_index);
    lexp.release();

    size_t bytes = translated->bytes();
    char *transexpand = reinterpret_cast<char *>(translated->release());

    // perform the set
    expansionlco_t lco{target, parms->n_digits};
    lco.contribute(bytes, transexpand);

    delete [] transexpand;
    return HPX_SUCCESS;
  }

  static void l_to_l_out_edge(Header *head, hpx_addr_t target, int n_digits) {
    // Get a parcel
    size_t msg_size = head->expansion_size + sizeof(LtoLParams);
    hpx_parcel_t *parc = hpx_parcel_acquire(nullptr, msg_size);
    assert(parc != nullptr);

    // Setup the parcel
    hpx_parcel_set_action(parc, l_to_l_);
    hpx_parcel_set_target(parc, target);

    // Fill the parcel payload
    LtoLParams *message = static_cast<LtoLParams *>(hpx_parcel_get_data(parc));
    message->bytes = head->expansion_size;
    message->n_digits = n_digits;
    memcpy(message->payload, head->payload, head->expansion_size);

    // Send the parcel - this will release the acquired parcel
    hpx_parcel_send_sync(parc);
  }

  static int l_to_l_handler(LtoLParams *parms, size_t UNUSED) {
    hpx_addr_t target = hpx_thread_current_target();

    void *workaround{nullptr};
    assert(hpx_gas_try_pin(target, &workaround));
    Header *trg_data =
        static_cast<Header *>(hpx_lco_user_get_user_data(workaround));
    Index t_index = trg_data->index;
    int to_child = t_index.which_child();
    LocalData<DomainGeometry> geo = trg_data->domain.value();
    hpx_gas_unpin(target);
    double t_size = geo->size_from_level(t_index.level());

    expansion_t lexp{parms->payload, parms->bytes, parms->n_digits};
    auto translated = lexp.L_to_L(to_child, t_size);
    lexp.release();

    size_t bytes = translated->bytes();
    char *transexpand = reinterpret_cast<char *>(translated->release());

    expansionlco_t total{target, parms->n_digits};
    total.contribute(bytes, transexpand);

    delete [] transexpand;

    return HPX_SUCCESS;
  }

  static void m_to_t_out_edge(Header *head, hpx_addr_t target, int n_digits) {
    LocalData<DomainGeometry> geo = head->domain.value();
    double scale = geo->size_from_level(head->index.level());

    // NOTE: we do not put in the correct number of targets. This is fine
    // because contribute_M_to_T does not rely on this information.
    targetlco_t destination{target, 0};
    destination.contribute_M_to_T(head->expansion_size, head->payload,
                                  n_digits, scale);
  }

  static void l_to_t_out_edge(Header *head, hpx_addr_t target, int n_digits) {
    LocalData<DomainGeometry> geo = head->domain.value();
    double scale = geo->size_from_level(head->index.level());

    // NOTE: we do not put in the correct number of targets. This is fine
    // because contribute_L_to_T does not rely on this information.
    targetlco_t destination{target, 0};
    destination.contribute_L_to_T(head->expansion_size, head->payload,
                                  n_digits, scale);
  }

  static int s_to_t_handler(Source *sources, int n_sources,
                            hpx_addr_t target, int n_trg) {
    targetlco_t targets{target, n_trg};
    targets.contribute_S_to_T(n_sources, sources);
    return HPX_SUCCESS;
  }

  static int create_from_expansion_handler(void *payload, size_t bytes) {
    hpx_addr_t gdata = hpx_lco_user_new(bytes, init_, operation_,
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
  static hpx_action_t s_to_m_;
  static hpx_action_t s_to_l_;
  static hpx_action_t m_to_l_;
  static hpx_action_t l_to_l_;
  static hpx_action_t s_to_t_;
  static hpx_action_t create_from_expansion_;

  hpx_addr_t data_;     // this is the LCO
  int n_digits_;
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
hpx_action_t ExpansionLCO<S, T, E, M, D>::s_to_m_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t ExpansionLCO<S, T, E, M, D>::s_to_l_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t ExpansionLCO<S, T, E, M, D>::m_to_l_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t ExpansionLCO<S, T, E, M, D>::l_to_l_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t ExpansionLCO<S, T, E, M, D>::s_to_t_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t ExpansionLCO<S, T, E, M, D>::create_from_expansion_
    = HPX_ACTION_NULL;


} // namespace dashmm


#endif // __DASHMM_EXPANSION_REF_H__
