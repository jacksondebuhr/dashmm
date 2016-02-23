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


/// \file include/targetlco.h
/// \brief TargetLCO object definition


#include <cstring>

#include <hpx/hpx.h>

#include "dashmm/targetref.h"


namespace dashmm {


/// Forward declaration of Evaluator so that we can become friends
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class> class Method>
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
/// This is a template class parameterized by the Source, Target, Expansion
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

  using targetref_t = TargetRef<Target>;

  /// Construct a default object
  TargetLCO() : lco_{HPX_NULL}, n_targs_{0} { }

  /// Construct from an existing LCO
  TargetLCO(hpx_addr_t data, int n_targs) : lco_{data}, n_targs_{n_targs} { }

  /// Construct an LCO from input TargetRef. This will create the LCO.
  explicit TargetLCO(const targetref_t &targets) {
    Init init{targets};
    lco_ = hpx_lco_user_new(sizeof(Data), init_, operation_, predicate_,
                            &init, sizeof(init));
    assert(lco_ != HPX_NULL);
    n_targs_ = targets.n();
  }

  /// Destroy the LCO. Use carefully!
  void destroy() {
    if (lco_ != HPX_NULL) {
      hpx_lco_delete_sync(lco_);
      lco_ = HPX_NULL;
    }
  }

  /// The global address of the referred to object
  hpx_addr_t lco() const {return lco_;}

  int n() const {return n_targs_;}

  /// Indicate to the underlying LCO that all operations have been scheduled
  void finalize() const {
    if (lco_ != HPX_NULL) {
      int code = kFinish;
      hpx_lco_set_lsync(lco_, sizeof(code), &code, HPX_NULL);
    }
  }

  /// Indicate to the underlying LCO that it should expect \param num more
  /// operations
  void schedule(int num) const {
    if (lco_ != HPX_NULL) {
      int input[2] = {kSetOnly, num};
      // NOTE: This one must be remote complete so that there is no timing issue
      // between scheduling an input and finalizing the schedule.
      hpx_lco_set_rsync(lco_, sizeof(int) * 2, input);
    }
  }

  /// Contribute a S->T operation to the referred targets
  ///
  /// \param n - the number of sources
  /// \param sources - the sources themselves
  void contribute_S_to_T(int n, source_t *sources) const {
    // NOTE: we assume this is called local to the sources
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
  /// \param scale - scaling factor
  void contribute_M_to_T(size_t bytes, void *data,
                         int n_digits, double scale) const {
    size_t inputsize = sizeof(MtoT) + bytes;
    MtoT *input = reinterpret_cast<MtoT *>(new char [inputsize]);
    assert(input);
    input->code = kMtoT;
    input->n_digits = n_digits;
    input->scale = scale;
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
  /// \param scale - scaling factor
  void contribute_L_to_T(size_t bytes, void *data,
                         int n_digits, double scale) const {
    size_t inputsize = sizeof(LtoT) + bytes;
    LtoT *input = reinterpret_cast<LtoT *>(new char [inputsize]);
    assert(input);
    input->code = kLtoT;
    input->n_digits = n_digits;
    input->scale = scale;
    input->bytes = bytes;
    memcpy(input->data, data, bytes);
    hpx_lco_set_lsync(lco_, inputsize, input, HPX_NULL);
  }

 private:
  /// Make Evaluator a friend -- it handles action registration
  friend class Evaluator<Source, Target, Expansion, Method>;

  /// LCO data type
  struct Data {
    int arrived;
    int scheduled;
    int finished;
    targetref_t targets;
  };

  /// LCO initialization type
  struct Init {
    targetref_t targets;
  };

  /// S->T parameters type
  struct StoT {
    int code;
    int count;
    source_t sources[];
  };

  /// M->T parameters type
  struct MtoT {
    int code;
    int n_digits;
    double scale;
    size_t bytes;
    char data[];
  };

  /// L->T parameters type
  struct LtoT {
    int code;
    int n_digits;
    double scale;
    size_t bytes;
    char data[];
  };

  /// Codes to define what the LCO set is doing.
  enum SetCodes {
    kSetOnly = 1,
    kStoT = 2,
    kMtoT = 3,
    kLtoT = 4,
    kFinish = 5
  };

  /// Initialize the LCO
  static void init_handler(Data *i, size_t bytes,
                           Init *init, size_t init_bytes) {
    i->arrived = 0;
    i->scheduled = 0;
    i->finished = 0;
    i->targets = init->targets;
  }

  /// The 'set' operation on the LCO
  ///
  /// This takes a number of forms based on the input code.
  static void operation_handler(Data *lhs, void *rhs, size_t bytes) {
    int *code = static_cast<int *>(rhs);
    if (*code == kSetOnly) {    // this is a pair of ints, a code and a count
      lhs->scheduled += code[1];
    } else if (*code == kStoT) {
      // The input contains the sources
      StoT *input = static_cast<StoT *>(rhs);

      // The LCO data contains the reference to the targets, which must be
      // pinned.
      target_t *targets{nullptr};
      // NOTE: This should succeed because we place the LCO at the same locality
      // as the actual target data.
      assert(hpx_gas_try_pin(lhs->targets.data(), (void **)&targets));

      expansion_t expand{nullptr, 0, -1};
      expand.S_to_T(input->sources, &input->sources[input->count],
                     targets, &targets[lhs->targets.n()]);
      expand.release(); // NOTE: This is not strictly needed

      hpx_gas_unpin(lhs->targets.data());

      lhs->arrived += 1;
    } else if (*code == kMtoT) {
      MtoT *input = static_cast<MtoT *>(rhs);

      // The LCO data contains the reference to the targets, which must be
      // pinned.
      target_t *targets{nullptr};
      // NOTE: This should succeed because we place the LCO at the same locality
      // as the actual target data.
      assert(hpx_gas_try_pin(lhs->targets.data(), (void **)&targets));

      expansion_t expand{input->data, input->bytes, input->n_digits};
      expand.M_to_T(targets, &targets[lhs->targets.n()], input->scale);
      expand.release();

      hpx_gas_unpin(lhs->targets.data());

      lhs->arrived += 1;
    } else if (*code == kLtoT) {
      LtoT *input = static_cast<LtoT *>(rhs);

      // The LCO data contains the reference to the targets, which must be
      // pinned.
      target_t *targets{nullptr};
      // NOTE: This should succeed because we place the LCO at the same locality
      // as the actual target data.
      assert(hpx_gas_try_pin(lhs->targets.data(), (void **)&targets));

      expansion_t expand{input->data, input->bytes, input->n_digits};
      expand.L_to_T(targets, &targets[lhs->targets.n()], input->scale);
      expand.release();

      hpx_gas_unpin(lhs->targets.data());

      lhs->arrived += 1;
    } else if (*code == kFinish) {
      lhs->finished = 1;
    } else {
      assert(0 && "Incorrect code to TargetLCO");
    }
  }

  /// The LCO is set if it has been finalized, and all scheduled operations
  /// have taken place.
  static bool predicate_handler(Data *i, size_t bytes) {
    return i->finished && (i->arrived == i->scheduled);
  }

  /// The global address of the LCO
  hpx_addr_t lco_;
  /// The number of targets represented by the LCO.
  int n_targs_;

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
