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


#ifndef __DASHMM_TARGET_LCO_H__
#define __DASHMM_TARGET_LCO_H__


/// \file include/dashmm/targetlco.h
/// \brief TargetLCO object definition


#include <cstring>

#include <hpx/hpx.h>

#include "dashmm/arrayref.h"
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


/// Target LCO
///
/// This LCO manages the concurrent contribution to the target data. In
/// principle, the data might be modified by the LCO and some other means
/// simultaneously, as the LCO cannot exclude on data exterior to the LCO
/// data itself. However, discipline inside DASHMM, and a mild discipline
/// by the user will obviate any worry.
///
/// The LCO contains a TargetRef for the targets it 'owns' even though it
/// cannot really be said to own the data. DASHMM will allocate this LCO
/// at the same locality as the data referred to by the TargetRef, so they
/// should always be co-located. Generally speaking, the user does not
/// interact with the object very often. Mostly they will pass objects of this
/// type to ExpansionLCO objects.
///
/// This is a template class parameterized by the Source, Target, Expansion,
/// Method and DistroPolicy types for a particular evaluation of DASHMM.
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class Method,
          typename DistroPolicy>
class TargetLCO {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Method<Source, Target, Expansion, DistroPolicy>;

  using targetref_t = ArrayRef<Target>;

  /// Construct a default object
  TargetLCO() : lco_{HPX_NULL}, n_targs_{0} { }

  /// Construct from an existing LCO
  TargetLCO(hpx_addr_t data, size_t n_targs) : lco_{data}, n_targs_{n_targs} { }

  /// Construct an LCO from input TargetRef. This will create the LCO.
  TargetLCO(size_t n_inputs, const targetref_t &targets, hpx_addr_t where) {
    Data init{static_cast<int>(n_inputs), targets};
    hpx_call_sync(where, create_, &lco_, sizeof(lco_), &init, sizeof(init));
    assert(lco_ != HPX_NULL);
    n_targs_ = targets.n();
  }

  /// Destroy the LCO. Use carefully!
  void destroy() {
    if (lco_ != HPX_NULL) {
      hpx_lco_delete_sync(lco_);
      lco_ = HPX_NULL;
      n_targs_ = 0;
    }
  }

  /// The global address of the referred to object
  hpx_addr_t lco() const {return lco_;}

  /// The number of targets in the referred object
  size_t n() const {return n_targs_;}

  /// Contribute a S->T operation to the referred targets
  ///
  /// \param n - the number of sources
  /// \param sources - the sources themselves
  void contribute_S_to_T(size_t n, source_t *sources) const {
    size_t inputsize = sizeof(StoT) + sizeof(source_t) * n;
    StoT *input = reinterpret_cast<StoT *>(new char [inputsize]);
    assert(input);
    input->code = kStoT;
    input->count = n;
    memcpy(input->sources, sources, sizeof(source_t) * n);

    hpx_lco_set_lsync(lco_, inputsize, input, HPX_NULL);

    delete [] input;
  }

  /// Contribute a M->T operation to the referred targets
  ///
  /// \param bytes - the size of the serialized expansion data
  /// \param data - the serialized expansion data
  /// \param n_digits - accuracy of the expansion
  void contribute_M_to_T(size_t bytes, void *data,
                         int n_digits) const {
    size_t inputsize = sizeof(MtoT) + bytes;
    MtoT *input = reinterpret_cast<MtoT *>(new char [inputsize]);
    assert(input);
    input->code = kMtoT;
    input->n_digits = n_digits;
    input->bytes = bytes;
    memcpy(input->data, data, bytes);
    hpx_lco_set_lsync(lco_, inputsize, input, HPX_NULL);
    delete [] input;
  }

  /// Contribute a L->T operation to the referred targets
  ///
  /// \param bytes - the size of the serialized expansion data
  /// \param data - the serialized expansion data
  /// \param n_digits - accuracy of the expansion
  void contribute_L_to_T(size_t bytes, void *data,
                         int n_digits) const {
    size_t inputsize = sizeof(LtoT) + bytes;
    LtoT *input = reinterpret_cast<LtoT *>(new char [inputsize]);
    assert(input);
    input->code = kLtoT;
    input->n_digits = n_digits;
    input->bytes = bytes;
    memcpy(input->data, data, bytes);
    hpx_lco_set_lsync(lco_, inputsize, input, HPX_NULL);
    delete [] input;
  }

 private:
  /// Make Evaluator a friend -- it handles action registration
  friend class Evaluator<Source, Target, Expansion, Method, DistroPolicy>;

  /// LCO data type
  struct Data {
    int yet_to_arrive;
    targetref_t targets;
  };

  /// S->T parameters type
  struct StoT {
    int code;
    size_t count;
    source_t sources[];
  };

  /// M->T parameters type
  struct MtoT {
    int code;
    int n_digits;
    size_t bytes;
    char data[];
  };

  /// L->T parameters type
  struct LtoT {
    int code;
    int n_digits;
    size_t bytes;
    char data[];
  };

  /// Codes to define what the LCO set is doing.
  enum SetCodes {
    kStoT = 0,
    kMtoT = 1,
    kLtoT = 2,
  };

  /// Initialize the LCO
  static void init_handler(Data *i, size_t bytes,
                           Data *init, size_t init_bytes) {
    *i = *init;
  }

  /// The 'set' operation on the LCO
  ///
  /// This takes a number of forms based on the input code.
  static void operation_handler(Data *lhs, void *rhs, size_t bytes) {
    int *code = reinterpret_cast<int *>(rhs);
    ReadBuffer readbuf{(char *)rhs, bytes};

    lhs->yet_to_arrive -= 1;
    assert(lhs->yet_to_arrive >= 0);

    if (*code == kStoT) {
      // The input contains the sources
      StoT *input = static_cast<StoT *>(rhs);

      // The LCO data contains the reference to the targets, which must be
      // pinned.
      target_t *targets{nullptr};
      // NOTE: This should succeed because we place the LCO at the same locality
      // as the actual target data.
      assert(hpx_gas_try_pin(lhs->targets.data(), (void **)&targets));

      expansion_t expand(ViewSet{});
      expand.S_to_T(input->sources, &input->sources[input->count],
                     targets, &targets[lhs->targets.n()]);

      hpx_gas_unpin(lhs->targets.data());
    } else if (*code == kMtoT) {
      // TODO: take a look at this MtoT type. We use nothing but the code here. 
      readbuf.interpret<MtoT>();

      // The LCO data contains the reference to the targets, which must be
      // pinned.
      target_t *targets{nullptr};
      // NOTE: This should succeed because we place the LCO at the same locality
      // as the actual target data.
      assert(hpx_gas_try_pin(lhs->targets.data(), (void **)&targets));

      ViewSet views{};
      views.interpret(readbuf);
      expansion_t expand(views);
      expand.M_to_T(targets, &targets[lhs->targets.n()]);
      expand.release();

      hpx_gas_unpin(lhs->targets.data());
    } else if (*code == kLtoT) {
      // TODO: take a look at this LtoT type. We use nothing but the code here. 
      readbuf.interpret<LtoT>(); 

      // The LCO data contains the reference to the targets, which must be
      // pinned.
      target_t *targets{nullptr};
      // NOTE: This should succeed because we place the LCO at the same locality
      // as the actual target data.
      assert(hpx_gas_try_pin(lhs->targets.data(), (void **)&targets));

      ViewSet views{};
      views.interpret(readbuf);
      expansion_t expand(views);
      expand.L_to_T(targets, &targets[lhs->targets.n()]);
      expand.release();

      hpx_gas_unpin(lhs->targets.data());
    } else {
      assert(0 && "Incorrect code to TargetLCO");
    }
  }

  /// The LCO is set if it has been finalized, and all scheduled operations
  /// have taken place.
  static bool predicate_handler(Data *i, size_t bytes) {
    return (i->yet_to_arrive == 0);
  }

  /// Action to create the LCO at a given locality
  static int create_at_locality(Data *init, size_t bytes) {
    hpx_addr_t retval = hpx_lco_user_new(bytes, init_, operation_, predicate_,
                                         init, bytes);
    return HPX_THREAD_CONTINUE(retval);
  }

  /// The global address of the LCO
  hpx_addr_t lco_;
  /// The number of targets represented by the LCO.
  size_t n_targs_;

  /// HPX function for LCO initialization
  static hpx_action_t init_;
  /// HPX function for LCO set
  static hpx_action_t operation_;
  /// HPX function for LCO predicate
  static hpx_action_t predicate_;
  /// Action for locality aware construction
  static hpx_action_t create_;
};


template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t TargetLCO<S, T, E, M, D>::init_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t TargetLCO<S, T, E, M, D>::operation_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t TargetLCO<S, T, E, M, D>::predicate_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t TargetLCO<S, T, E, M, D>::create_ = HPX_ACTION_NULL;


} // namespace dashmm


#endif // __DASHMM_TARGET_LCO_H__
