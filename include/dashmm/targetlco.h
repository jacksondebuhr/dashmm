// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2017, Trustees of Indiana University,
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


/// \file
/// \brief TargetLCO object definition


#include <cstring>

#include <hpx/hpx.h>

#include "dashmm/arrayref.h"
#include "dashmm/traceevents.h"
#include "dashmm/viewset.h"


namespace dashmm {


/// Forward declaration of TargetLCORegistrar so that we can become friends
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class> class Method>
class TargetLCORegistrar;


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
/// and Method types for a particular evaluation of DASHMM.
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class> class Method>
class TargetLCO {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Method<Source, Target, Expansion>;

  using targetref_t = ArrayRef<Target>;

  /// Construct a default object
  TargetLCO() : lco_{HPX_NULL} { }

  /// Construct from an existing LCO
  TargetLCO(hpx_addr_t data) : lco_{data} { }

  /// Construct an LCO from input TargetRef
  ///
  /// This will construct the LCO at the locality specified by @p where.
  /// The LCO will expect @p n_inputs contributions.
  ///
  /// \param n_inputs - the number of inputs to the LCO
  /// \param targets - ArrayRef indicating the global memory that the LCO is
  ///                  representing
  TargetLCO(size_t n_inputs, const targetref_t &targets) {
    Data init{static_cast<int>(n_inputs), targets};
    lco_ = hpx_lco_user_new(sizeof(init), init_, operation_,
                            predicate_, &init, sizeof(init));
    assert(lco_ != HPX_NULL);
  }

  /// Destroy the LCO
  void destroy() {
    if (lco_ != HPX_NULL) {
      hpx_lco_delete_sync(lco_);
      lco_ = HPX_NULL;
    }
  }

  /// The global address of the referred to object
  hpx_addr_t lco() const {return lco_;}

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
    if (n != 0) {
      memcpy(input->sources, sources, sizeof(source_t) * n);
    }

    hpx_lco_set_lsync(lco_, inputsize, input, HPX_NULL);

    delete [] input;
  }

  /// Contribute a M->T operation to the referred targets
  ///
  /// \param bytes - the size of the serialized expansion data
  /// \param data - the serialized expansion data
  void contribute_M_to_T(size_t bytes, void *data) const {
    size_t inputsize = sizeof(MtoT) + bytes;
    MtoT *input = reinterpret_cast<MtoT *>(new char [inputsize]);
    assert(input);
    input->code = kMtoT;
    input->bytes = bytes;
    memcpy(input->data, data, bytes);
    hpx_lco_set_lsync(lco_, inputsize, input, HPX_NULL);
    delete [] input;
  }

  /// Contribute a L->T operation to the referred targets
  ///
  /// \param bytes - the size of the serialized expansion data
  /// \param data - the serialized expansion data
  void contribute_L_to_T(size_t bytes, void *data) const {
    size_t inputsize = sizeof(LtoT) + bytes;
    LtoT *input = reinterpret_cast<LtoT *>(new char [inputsize]);
    assert(input);
    input->code = kLtoT;
    input->bytes = bytes;
    memcpy(input->data, data, bytes);
    hpx_lco_set_lsync(lco_, inputsize, input, HPX_NULL);
    delete [] input;
  }

 private:
  /// Become friends with TargetLCORegistrar
  friend class TargetLCORegistrar<Source, Target, Expansion, Method>;

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
    size_t bytes;
    char data[];
  };

  /// L->T parameters type
  struct LtoT {
    int code;
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

    if (lhs->targets.data() == nullptr) {
      return;
    }

    if (*code == kStoT) {
      // The input contains the sources
      EVENT_TRACE_DASHMM_STOT_BEGIN();
      StoT *input = static_cast<StoT *>(rhs);

      if (input->count) {
        target_t *targets{lhs->targets.data()};

        expansion_t expand(ViewSet{});
        expand.S_to_T(input->sources, &input->sources[input->count],
                       targets, &targets[lhs->targets.n()]);

      }
      EVENT_TRACE_DASHMM_STOT_END();
    } else if (*code == kMtoT) {
      EVENT_TRACE_DASHMM_MTOT_BEGIN();
      readbuf.interpret<MtoT>();

      target_t *targets{lhs->targets.data()};

      ViewSet views{};
      views.interpret(readbuf);
      expansion_t expand(views);
      expand.M_to_T(targets, &targets[lhs->targets.n()]);
      expand.release();
      EVENT_TRACE_DASHMM_MTOT_END();
    } else if (*code == kLtoT) {
      EVENT_TRACE_DASHMM_LTOT_BEGIN();
      readbuf.interpret<LtoT>();

      target_t *targets{lhs->targets.data()};

      ViewSet views{};
      views.interpret(readbuf);
      expansion_t expand(views);
      expand.L_to_T(targets, &targets[lhs->targets.n()]);
      expand.release();
      EVENT_TRACE_DASHMM_LTOT_END();
    } else {
      assert(0 && "Incorrect code to TargetLCO");
    }
  }

  /// The LCO is set if all scheduled operations have taken place.
  static bool predicate_handler(Data *i, size_t bytes) {
    return (i->yet_to_arrive == 0);
  }

  /// The global address of the LCO
  hpx_addr_t lco_;

  /// HPX function for LCO initialization
  static hpx_action_t init_;
  /// HPX function for LCO set
  static hpx_action_t operation_;
  /// HPX function for LCO predicate
  static hpx_action_t predicate_;
};


template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t TargetLCO<S, T, E, M>::init_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t TargetLCO<S, T, E, M>::operation_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t TargetLCO<S, T, E, M>::predicate_ = HPX_ACTION_NULL;


} // namespace dashmm


#endif // __DASHMM_TARGET_LCO_H__
