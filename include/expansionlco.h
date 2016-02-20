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


#ifndef __DASHMM_EXPANSION_REF_H__
#define __DASHMM_EXPANSION_REF_H__


/// \file include/expansionlco.h
/// \brief Interface to Expansion LCO


#include <memory>
#include <vector>

#include <hpx/hpx.h>

#include "include/index.h"
#include "include/point.h"
#include "include/sourceref.h"
#include "include/targetlco.h"
#include "include/targetref.h"
#include "include/types.h"


namespace dashmm {


/// Forward declaration of Evaluator so that we can become friends
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename, typename> class Method>
class Evaluator<Source, Target, Expansion, Method>;


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
          template <typename, typename, typename> class Method>
class ExpansionLCO {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Method<Source, Target, expansion_t>;

  using sourceref_t = SourceRef<Source>;
  using targetref_t = TargetRef<Target>;
  using targetlco_t = TargetLCO<Source, Target, Expansion, Method>;

  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, Method>;

  /// Construct the expansion from a given global address.
  ExpansionLCO(hpx_addr_t addr, int n_digits)
      : data_{addr}, n_digits_{n_digits} { }

  /// Construct from an existing expansion - this will create a new LCO
  ExpansionLCO(std::unique_ptr<expanstion_t> exp, hpx_addr_t where) {
    if (exp == nullptr) {
      return ExpansionLCO{HPX_NULL, -1};
    }

    size_t bytes = exp->bytes();
    int n_digits = exp->accuracy();
    char *ldata = reinterpret_cast<char *>(exp->release());

    hpx_addr_t retval{HPX_NULL};
    hpx_call_sync(where, create_from_expansion_, &retval, sizeof(retval),
                  ldata, bytes);
    delete [] ldata;

    return ExpansionLCO{retval, n_digits};
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
  /// Note that this sets the expansion to be equal to the computed
  /// multipole expansion. Further, this must not be called
  /// after finalize().
  ///
  /// \param center - the center of the computed expansion
  /// \param sources - a reference to the sources from which to compute the
  ///                  multipole expansion
  /// \param scale - scaling factor
  void S_to_M(Point center, sourceref_t sources, double scale) const {
    schedule();
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
  /// current contents of the expansion. Further, this must not be called
  /// after finalize().
  ///
  /// \param center - the center of the computed expansion
  /// \param sources - a reference to the sources from which to compute the
  ///                  local expansion
  /// \param scale - scaling factor
  void S_to_L(Point center, sourceref_t sources, double scale) const {
    schedule();
    int nsrc = sources.n();
    double cx = center.x();
    double cy = center.y();
    double cz = center.z();
    hpx_call(sources.data(), s_to_l_, HPX_NULL, &nsrc, &cx, &cy, &cz, &scale,
             &data_, &n_digits_);
  }

  /// Contribute a translated multipole moment to this expansion
  ///
  /// This will translate the given multipole expansion and then add it to
  /// this expansion. This is an asynchronous operation; this can return
  /// before the contribution to this expansion has been made. This must not
  /// be called after finalize().
  ///
  /// \param source - the multipole expansion to translate
  /// \param from_child - the child from which the expansion will occur
  /// \param s_size - the size of the child node for @p source
  void M_to_M(expansionlco_t source, int from_child, double s_size) const {
    schedule();
    hpx_call_when(source.data(), source.data(), m_to_m_, HPX_NULL,
                  &data_, &n_digits_, &from_child, &s_size);
  }

  /// Contribute a translated multipole moment to this local expansion
  ///
  /// This will translate a multipole expansion into a local expansion and
  /// add it to this expansion. This operation is asynchronous; this can
  /// return before the contribution has been made. This must not be called
  /// after finalize().
  ///
  /// \param source - the multipole expansion to translate
  /// \param s_index - the index of the node containing @p source
  /// \param s_size - the size of the node containing @p source
  /// \param t_index - the index of the node containing this expansion
  void M_to_L(expansionlco_t source, Index s_index, double s_size,
              Index t_index) const {
    schedule();
    MtoLParams args{*this, s_index, s_size, t_index, n_digits_};
    hpx_call_when(source.data(), source.data(), m_to_l_, HPX_NULL,
                  &args, sizeof(args));
  }

  /// Contribute a translated local expansion to this local expansion
  ///
  /// This will translate the given local expansion into a local expansion
  /// that can be added to this expansion. This operation is asynchronous;
  /// this can return before the contribution to this expansion has been made.
  /// This must not be called after finalize().
  ///
  /// \param source - the local expansion to translate
  /// \param to_child - the child to which the expansion is being translated
  /// \param t_size - the size of the child node
  void L_to_L(expansionlco_t source, int to_child, double t_size) const {
    schedule();
    hpx_call_when(source.data(), source.data(), l_to_l_, HPX_NULL,
                  &data_, &n_digits_, &to_child, &t_size);
  }

  /// Apply the effect of a multipole expansion to targets
  ///
  /// This will compute the effect of this multipole expansion on the given
  /// @p targets. Note that this is an asynchronous operation. This will only
  /// perform work once this expansion is ready.
  ///
  /// \param targets - the target for which the multipole expansion is applied
  /// \param scale - scaling factor
  void M_to_T(targetlco_t targets, double scale) const {
    targets.schedule(1);
    hpx_addr_t tsend = targets.data();
    hpx_call_when(data_, data_, m_to_t_, HPX_NULL, &n_digits_, &scale, &tsend);
  }

  /// Apply the effect of a local expansion to targets
  ///
  /// This will compute the effect of this local expansion on the given
  /// @p targets. Note that this is an asynchronous operation. This will
  /// only perform work once this expansion is ready.
  ///
  /// \param targets - the target for which the local expansion is applied
  /// \param scale - scaling factor
  void L_to_T(targetlco_t targets, double scale) const {
    targets.schedule(1);
    hpx_addr_t tsend = targets.data();
    hpx_call_when(data_, data_, l_to_t_, HPX_NULL, &n_digits_, &scale, &tsend);
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
    targets.schedule(1);
    int n_src = sources.n();
    hpx_addr_t tsend = targets.data();
    hpx_call(sources.data(), s_to_t_, HPX_NULL, &n_src, &tsend);
  }

  /// Add the given expansion to this expansion
  ///
  /// This will add the @p summand to this expansion. This is an asynchronous
  /// operation, and will only complete once @p summand is set. This does
  /// schedule a contribution, so this should only be called before the call
  /// to finalize(). This routine takes ownership of the supplied expansion,
  /// and will free any resources associated with the expansion.
  ///
  /// \param summand - a reference to the expansion to add to this one
  void add_expansion(expansionlco_t summand) {
    schedule();  // we are going to have another contribution
    hpx_call_when(summand.data(), summand.data(), add_,
                  HPX_NULL, &data_, &n_digits_);
  }

  /// Create a new expansion of the same type referred to by this object.
  ///
  /// \param center - the center point for the next expansion.
  /// \param n_digits - the accuracy parameter for the expansion
  ///
  /// \returns - the resulting expansion.
  std::unique_ptr<expansion_t> get_new_expansion(Point center,
                                                 int n_digits) const {
    return std::unique_ptr<expansion_t>{new expansion_t{center, n_digits}};
  }

  /// Signal to the expansion that all operations have been scheduled.
  ///
  /// This should be called only after all possible contributions to the
  /// expansion have been scheduled. The underlying LCO cannot trigger until
  /// finalize() has been called, and so any work dependent on this expansion
  /// would be on-hold until the call to finalize().
  void finalize() const {
    if (data_ != HPX_NULL) {
      int code = kFinish;
      hpx_lco_set_lsync(data_, sizeof(code), &code, HPX_NULL);
    }
  }

  /// Signal to the expansion that it should expect an operation
  ///
  /// This will inform the underlying LCO of another eventual contribution to
  /// its value. The asynchronous nature of the computation in DASHMM means that
  /// contributions will happen when they are ready. schedule() will inform
  /// that a contribution is on the way.
  void schedule() const {
    if (data_ != HPX_NULL) {
      int code = kSchedule;
      hpx_lco_set_rsync(data_, sizeof(code), &code);
    }
  }


  /// Contribute to the referred expansion
  ///
  /// This will setup the given @p payload with the correct internal code
  /// and will call the appropriate set operation on the referred LCO. This
  /// will result in the add_expansion method of the expansion being called.
  ///
  /// \param bytes - the size of the input serialized expansion
  /// \param payload - the serialized expansion data
  void contribute(size_t bytes, char *payload) {
    int *code = reinterpret_cast<int *>(payload);
    *code = kContribute;
    hpx_lco_set_lsync(data_, bytes, payload, HPX_NULL);
  }


 private:
  // Give the evaluator access so that it might register our actions
  friend class Evaluator<Source, Target, Expansion, Method>;

  ///////////////////////////////////////////////////////////////////
  // Types used internally
  ///////////////////////////////////////////////////////////////////

  /// Part of the internal representation of the Expansion LCO
  ///
  /// Expansion are user-defined LCOs. The data they contain are this object
  /// and the serialized expansion. This object gives the number of expected
  /// inputs, the number that have actually occurred, and a flag to indicate
  /// if all of the expected inputs have been scheduled.
  struct Header {
    int arrived;
    int scheduled;
    int finished;
    int payload_size;
    char payload[];
  };

  /// Behavior codes for the Expansion LCO
  ///
  /// The set operation for the Expansion LCO takes three forms. Two are simple:
  /// incrementing the number count of scheduled inputs, and finalizing the
  /// inputs. The third is for actually making contributions to the Expansion
  /// data.
  enum SetCodes {
    kFinish = 1,
    kSchedule = 2,
    kContribute = 3
  };

  // Marshalled parameter type for m_to_l_
  struct MtoLParams {
    expansionlco_t total;
    Index s_index;
    double s_size;
    Index t_index;
    int n_digits;
  };


  ///////////////////////////////////////////////////////////////////
  // LCO Implementation
  ///////////////////////////////////////////////////////////////////

  /// Initialization handler for Expansion LCOs
  ///
  /// This initialized an Expansion LCO given an input serialized
  /// expansion. Often, this will just be the default constructed expansion,
  /// but might be otherwise in specific cases.
  static void init_handler(Header *head, size_t bytes,
                           void *init, size_t init_bytes) {
    assert(bytes == init_bytes + sizeof(Header));
    head->arrived = 0;
    head->scheduled = 0;
    head->finished = 0;
    head->payload_size = init_bytes;
    memcpy(head->payload, init, init_bytes);
  }

  /// The set operation handler for the Expansion LCO
  ///
  /// This takes one of three forms. The input to this is either a single integer
  /// or a serialized Expansion. In the latter case, the reserved data at the
  /// beginning of the expansion serialization is used to give the operation
  /// code for the set.
  static void operation_handler(Header *lhs, void *rhs, size_t bytes) {
    int *code = static_cast<int *>(rhs);

    if (*code == kFinish) {
      assert(lhs->finished == 0);
      lhs->finished = 1;
    } else if (*code == kSchedule) {
      assert(lhs->finished == 0);
      lhs->scheduled += 1;
    } else if (*code == kContribute) {
      // increment the counter
      lhs->arrived += 1;
      if (lhs->finished) {
        assert(lhs->arrived <= lhs->scheduled);
      }

      // This is for the S->M contribution. That one occurs in place, and so
      // we need to contribute without actually doing anything.
      if (bytes <= sizeof(int)) {
        return;
      }

      // create an expansion from the rhs
      expansion_t::contents_t *input =
        static_cast<expansion_t::contents_t *>(rhs);
      // The start of contents_t must be two integers, the second being
      // n_digits.
      expansion_t incoming{input, bytes, code[1]};

      // create the expansion from the payload
      int *n_digits = reinterpret_cast<int *>(lhs->payload + sizeof(int));
      expansion_t expand{lhs->payload, lhs->payload_size, *n_digits};

      // add the one to the other
      expand.add_expansion(incoming);

      // release the data, because these objects do not actually own those buffers
      expand.release();
      incoming.release();

    } else {
      assert(0 && "Incorrect code to expansion LCO");
    }
  }

  /// The predicate to detect triggering of the Expansion LCO
  ///
  /// The expansion LCO is triggered if it has been finalized and the number
  /// of contributions match the number of scheduled operations.
  static bool predicate_handler(Header *i, size_t bytes) {
    return (i->finished && (i->arrived == i->scheduled));
  }


  ///////////////////////////////////////////////////////////////////
  // Other related actions
  ///////////////////////////////////////////////////////////////////

  static int s_to_m_handler(Source *sources, int n_src, double cx, double cy,
                            double cz, double scale, hpx_addr_t expand,
                            int n_digits) {
    void *temp{nullptr};

    // SMP assumption here: source and expansion LCO are on the same locality
    hpx_gas_try_pin(expand, &temp);
    Header *data = static_cast<Header *>(hpx_lco_user_get_user_data(temp));
    expansion_t local{data->payload, data->payload_size, n_digits};
    local->S_to_M(Point{cx, cy, cz}, sources, &sources[n_src], scale);
    local->release();
    hpx_gas_unpin(expand);

    int code = kContribute;
    hpx_lco_set_lsync(expand, sizeof(code), &code, HPX_NULL);

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

  static int m_to_m_handler(hpx_addr_t expand, int n_digits, int from_child,
                            double s_size) {
    hpx_addr_t target = hpx_thread_current_target();
    // HACK: This action is local to the expansion, so we getref here with
    // whatever as the size and things are okay...
    Header *ldata{nullptr};
    hpx_lco_getref(target, 1, (void **)&ldata);
    expansion_t lexp{ldata->payload, ldata->payload_size, n_digits};
    auto translated = lexp->M_to_M(from_child, s_size);
    lexp->release();
    hpx_lco_release(target, ldata);

    size_t bytes = translated->bytes();
    char *transexpand = reinterpret_cast<char *>(translated->release());

    expansionlco_t total{expand, n_digits};
    total.contribute(bytes, transexpand);

    delete [] transexpand;

    return HPX_SUCCESS;
  }

  static int m_to_l_handler(MtoLParams *parms, size_t UNUSED) {
    hpx_addr_t target = hpx_thread_current_target();
    // HACK: This action is local to the expansion, so we getref here with
    // whatever as the size and things are okay...
    Header *ldata{nullptr};
    hpx_lco_getref(target, 1, (void **)&ldata);

    expansion_t lexp{ldata->payload, ldata->payload_size, parms->n_digits};
    auto translated = lexp->M_to_L(parms->s_index, parms->s_size,
                                   parms->t_index);
    lexp->release();
    hpx_lco_release(target, ldata);

    size_t bytes = translated->bytes();
    char *transexpand = reinterpret_cast<char *>(translated->release());

    parms->total.contribute(bytes, transexpand);

    delete [] transexpand;

    return HPX_SUCCESS;
  }

  static int l_to_l_handler(hpx_addr_t expand, int n_digits, int to_child,
                            double t_size) {
    hpx_addr_t target = hpx_thread_current_target();
    // HACK: This action is local to the expansion, so we getref here with
    // whatever as the size and things are okay...
    Header *ldata{nullptr};
    hpx_lco_getref(target, 1, (void **)&ldata);
    expansion_t lexp{ldata->payload, ldata->payload_size, n_digits};
    auto translated = lexp->L_to_L(to_child, t_size);
    lexp->release();
    hpx_lco_release(target, ldata);

    size_t bytes = translated->bytes();
    char *transexpand = reinterpret_cast<char *>(translated->release());

    expansionlco_t total{expand, n_digits};
    total.contribute(bytes, transexpand);

    delete [] transexpand;

    return HPX_SUCCESS;
  }

  static int m_to_t_handler(int n_digits, double scale, hpx_addr_t targ) {
    targetlco_t targets{targ};
    // HACK: This action is local to the expansion, so we getref here with
    // whatever as the size and things are okay...
    Header *ldata{nullptr};
    hpx_lco_getref(hpx_thread_current_target(), 1, (void **)&ldata);
    targets.contribute_M_to_T(ldata->payload_size, ldata->payload,
                              n_digits, scale);
    hpx_lco_release(hpx_thread_current_target(), ldata);

    return HPX_SUCCESS;
  }

  static int l_to_t_handler(int n_digits, double scale, hpx_addr_t targ) {
    targetlco_t targets{targ};
    // HACK: This action is local to the expansion, so we getref here with
    // whatever as the size and things are okay...
    Header *ldata{nullptr};
    hpx_lco_getref(hpx_thread_current_target(), 1, (void **)&ldata);
    targets.contribute_L_to_T(ldata->payload_size, ldata->payload,
                              n_digits, scale);
    hpx_lco_release(hpx_thread_current_target(), ldata);

    return HPX_SUCCESS;
  }

  static int s_to_t_handler(Source *sources, int n_sources,
                            hpx_addr_t target) {
    targetlco_t targets{target};
    targets.contribute_S_to_T(n_sources, sources);
    return HPX_SUCCESS;
  }

  static int add_handler(hpx_addr_t expand, int n_digits) {
    expansionlco_t total{expand, n_digits};

    Header *ldata{nullptr};
    // HACK: This action is local to the expansion, so we getref here with
    // whatever as the size and things are okay...
    hpx_addr_t target = hpx_thread_current_target();
    hpx_lco_getref(target, 1, (void **)&ldata);
    total.contribute(ldata->payload_size, ldata->payload);
    hpx_lco_release(target, ldata);

    return HPX_SUCCESS;
  }

  static int create_from_expansion_handler(void *payload, size_t bytes) {
    size_t total_size = sizeof(Header) + bytes;
    hpx_addr_t gdata = hpx_lco_user_new(total_size, init_, operation_,
                                        predicate_, payload, bytes);
    assert(gdata != HPX_NULL);
    return HPX_THREAD_CONTINUE(gdata);
  }


  // The functions implementing the user LCO
  static hpx_action_t init_;
  static hpx_action_t operation_;
  static hpx_action_t predicate_;

  // The actions for the various operations
  static hpx_action_t s_to_m_;
  static hpx_action_t s_to_l_;
  static hpx_action_t m_to_m_;
  static hpx_action_t m_to_l_;
  static hpx_action_t l_to_l_;
  static hpx_action_t m_to_t_;
  static hpx_action_t l_to_t_;
  static hpx_action_t s_to_t_;
  static hpx_action_t add_;
  static hpx_action_t create_from_expansion_;

  hpx_addr_t data_;     // this is the LCO
  int n_digits_;
};

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename, typename> class M>
hpx_action_t ExpansionRef<S, T, E, M>::init_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename, typename> class M>
hpx_action_t ExpansionRef<S, T, E, M>::operation_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename, typename> class M>
hpx_action_t ExpansionRef<S, T, E, M>::predicate_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename, typename> class M>
hpx_action_t ExpansionRef<S, T, E, M>::s_to_m_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename, typename> class M>
hpx_action_t ExpansionRef<S, T, E, M>::s_to_l_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename, typename> class M>
hpx_action_t ExpansionRef<S, T, E, M>::m_to_m_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename, typename> class M>
hpx_action_t ExpansionRef<S, T, E, M>::m_to_l_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename, typename> class M>
hpx_action_t ExpansionRef<S, T, E, M>::l_to_l_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename, typename> class M>
hpx_action_t ExpansionRef<S, T, E, M>::m_to_t_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename, typename> class M>
hpx_action_t ExpansionRef<S, T, E, M>::l_to_t_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename, typename> class M>
hpx_action_t ExpansionRef<S, T, E, M>::s_to_t_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename, typename> class M>
hpx_action_t ExpansionRef<S, T, E, M>::add_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename, typename> class M>
hpx_action_t ExpansionRef<S, T, E, M>::create_from_expansion_ = HPX_ACTION_NULL;


} // namespace dashmm


#endif // __DASHMM_EXPANSION_REF_H__
