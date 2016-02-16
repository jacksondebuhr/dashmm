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


#ifndef __DASHMM_PARTICLE_H__
#define __DASHMM_PARTICLE_H__


/// \file include/particle.h
/// \brief Source and Target particle types.

#include <cassert>
#include <cstring>

#include <complex>

#include <hpx/hpx.h>

#include "include/point.h"
#include "include/types.h"


namespace dashmm {


/// Source concept in DASHMM
///
/// To qualify as a Source for DASHMM, a type should be trivially copyable,
/// and should provide at least two members accessible by name:
///
/// Point position;
/// double charge;
///
///
/// NOTE: For specific cases, you can add further data to that used by the
/// Expansions. The member will need to exist, or there will be some compilation
/// error.



/// Reference to a set of sources
///
/// This is a reference object, meaning that it refers to the Source data in
/// the GAS, but does not contain those data. As such, one can pass this
/// class by value without worry.
template <typename Source>
class SourceRef {
 public:
  using source_t = Source;

  /// Default constructor.
  SourceRef()
      : data_{HPX_NULL}, n_{0}, n_tot_{0} { }

  /// Construct from a specific address and count.
  SourceRef(hpx_addr_t data, size_t n, size_t n_tot)
      : data_{data}, n_{n}, n_tot_{n_tot} { }

  /// Destroy the particle data in GAS.
  ///
  /// This is needed because this object is a reference, and the destruction of
  /// this object only destroys the reference.
  ///
  /// Note that slices of another SourceRef should not be destroyed.
  /// TODO: Add some logic to prevent that? Another member is_slice_?
  void destroy() {
    if (data_ != HPX_NULL) {
      hpx_gas_free_sync(data_);
      data_ = HPX_NULL;
      n_ = 0;
    }
  }

  /// Returns if the reference is valid
  bool valid() const {return data_ != HPX_NULL;}

  /// Get a reference to a slide of the current reference
  ///
  /// This will return a SourceRef to a consecutive chunk of the records that
  /// begin at an offset from the start of this reference and that will contain
  /// n entries. If the input arguments to this method are invalid, then an
  /// invalid reference will be returned.
  ///
  /// \param offset - the offset from the start of this reference
  /// \param n - the number of entries in the slice
  ///
  /// \returns - the resulting SourceRef; may be invalid.
  SourceRef slice(size_t offset, size_t n) const {
    if (offset > n_) {
      return SourceRef{HPX_NULL, 0, 0};
    }
    if (offset + n > n_) {
      return SourceRef{HPX_NULL, 0, 0};
    }
    hpx_addr_t addr = hpx_addr_add(data_, sizeof(Source) * offset,
                                   sizeof(Source) * n_tot_);
    return SourceRef{addr, n, n_tot_};
  }

  /// Returns the number of Source records referred to.
  size_t n() const {return n_;}

  /// Returns the total number of Source records.
  size_t n_tot() const {return n_tot_;}

  /// Returns the global address of the referred to data.
  hpx_addr_t data() const {return data_;}

 private:
  hpx_addr_t data_;
  size_t n_;
  size_t n_tot_;
};



/// Target concept in DASHMM
///
/// To work as a Target type in DASHMM, a type must be trivially copyable and
/// must provide two members:
///
/// Point position
/// std::complex<double> phi
///
/// NOTE: Specific expansions may change these requirements. Please see the
/// documentation about the Expansion of interest.


/// Reference to a set of targets
///
/// This is a reference object, meaning that it refers to the Target data in
/// the GAS, but does not contain those data. As such, one can pass this
/// class by value without worry.
template <typename Target>
class TargetRef {
 public:
  using target_t = Target;

  /// Default constructor.
  TargetRef()
      : data_{HPX_NULL}, n_{0}, n_tot_{0} { }

  /// Construct from a specific address and count.
  TargetRef(hpx_addr_t data, size_t n, size_t n_tot)
      : data_{data}, n_{n}, n_tot_{n_tot} { }

  /// Destroy the particle data in GAS.
  ///
  /// This is needed because this object is a reference, and the destruction of
  /// this object only destroys the reference.
  ///
  /// Note that slices of another TargetRef should not be destroyed.
  /// TODO: Add some logic to prevent that? Another member is_slice_?
  void destroy() {
    if (data_ != HPX_NULL) {
      hpx_gas_free_sync(data_);
      data_ = HPX_NULL;
      n_ = 0;
    }
  }

  /// Returns if the reference is valid
  bool valid() const {return data_ != HPX_NULL;}

  /// Get a reference to a slide of the current reference
  ///
  /// This will return a TargetRef to a consecutive chunk of the records that
  /// begin at an offset from the start of this reference and that will contain
  /// n entries. If the input arguments to this method are invalid, then an
  /// invalid reference will be returned.
  ///
  /// \param offset - the offset from the start of this reference
  /// \param n - the number of entries in the slice
  ///
  /// \returns - the resulting TargetRef; may be invalid.
  TargetRef slice(size_t offset, size_t n) const {
    if (offset > n_) {
      return TargetRef{HPX_NULL, 0, 0};
    }
    if (offset + n > n_) {
      return TargetRef{HPX_NULL, 0, 0};
    }
    hpx_addr_t addr = hpx_addr_add(data_, sizeof(Target) * offset,
                                   sizeof(Target) * n_tot_);
    return TargetRef{addr, n, n_tot_};
  }

  /// Returns the number of Source records referred to.
  size_t n() const {return n_;}

  /// Returns the total number of Source records.
  size_t n_tot() const {return n_tot_;}

  /// Returns the global address of the referred to data.
  hpx_addr_t data() const {return data_;}

 private:
  hpx_addr_t data_;
  size_t n_;
  size_t n_tot_;
};


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
/// NOTE: We need all four here because otherwise, we would possibly have two
/// methods for the same Source, Target, Expansion triple try to register the
/// same functions. So that is no good. So we have to have all four.
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename, typename> class Method>
class TargetLCO {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Method<Source, Target, expansion_t>;

  using targetref_t = TargetRef<Target>;

  /// Construct a default object
  TargetLCO() : lco_{HPX_NULL} { }

  /// Construct from an existing LCO
  explicit TargetLCO(hpx_addr_t data) : lco_{data} { }

  /// Construct an LCO from input TargetRef. This will create the LCO.
  explicit TargetRef(targetref_t targets) {
    Init init{targets};
    lco_ = hpx_lco_user_new(sizeof(Data), init_, operation_, predicate_,
                            &init, sizeof(init));
    assert(lco_ != HPX_NULL);
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
  void contribute_S_to_T(int n, Source *sources) const {
    // NOTE: we assume this is called local to the sources
    size_t inputsize = sizeof(StoT) + sizeof(Source) * n;
    StoT *input = reinterpret_cast<StoT *>(new char [inputsize]);
    assert(input);
    input->code = kStoT;
    input->count = n;
    memcpy(input->sources, source, sizeof(Source) * n);

    hpx_lco_set_lsync(lco_, inputsize, input, HPX_NULL);

    delete [] input;
  }

  /// Contribute a M->T operation to the referred targets
  ///
  /// \param bytes - the size of the serialized expansion data
  /// \param data - the serialized expansion data
  /// \param n_digits - accuracy of the expansion
  /// \param scale - scaling factor
  void contribute_M_to_T(size_t bytes, expansion_t::contents_t *data,
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
  void contribute_L_to_T(size_t bytes, expansion_t::contents_t *data,
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
  // Make Evaluator a friend -- it handles action registration
  friend class Evaluator<Source, Target, Expansion, Method>;

  // some types used by the implementation
  struct Data {
    int arrived;
    int scheduled;
    int finished;
    targetref_t targets;
  };

  struct Init {
    targetref_t targets;
  };

  struct StoT {
    int code;
    int count;
    Source sources[];
  };

  struct MtoT {
    int code;
    int n_digits;
    double scale;
    size_t bytes;
    char data[];
  };

  struct LtoT {
    int code;
    int n_digits;
    double scale;
    size_t bytes;
    char data[];
  };

  enum SetCodes {
    kSetOnly = 1,
    kStoT = 2,
    kMtoT = 3,
    kLtoT = 4,
    kFinish = 5
  };

  // The functions implementing the LCO
  static void init_handler(Data *i, size_t bytes,
                           Init *init, size_t init_bytes) {
    i->arrived = 0;
    i->scheduled = 0;
    i->finished = 0;
    i->targets = init->targets;
  }

  static void operation_handler(Data *lhs, void *rhs, size_t bytes) {
    int *code = static_cast<int *>(rhs);
    if (*code == kSetOnly) {    // this is a pair of ints, a code and a count
      lhs->scheduled += code[1];
    } else if (*code == kStoT) {
      // The input contains the sources
      StoT *input = static_cast<StoT *>(rhs);

      // The LCO data contains the reference to the targets, which must be
      // pinned.
      target_t *targets{nullptr}
      // NOTE: This should succeed because we place the LCO at the same locality
      // as the actual target data.
      assert(hpx_gas_try_pin(lhs->targets.data(), (void **)targets));

      expansion_t expand{nullptr, 0, -1};
      expand->S_to_T(input->sources, &input->sources[input->count],
                     targets, &targets[lhs->targets.n()]);
      expansion->release(); // NOTE: This is not strictly needed

      hpx_gas_unpin(lhs->targets.data());

      lhs->arrived += 1;
    } else if (*code == kMtoT) {
      MtoT *input = static_cast<MtoT *>(rhs);

      // The LCO data contains the reference to the targets, which must be
      // pinned.
      target_t *targets{nullptr}
      // NOTE: This should succeed because we place the LCO at the same locality
      // as the actual target data.
      assert(hpx_gas_try_pin(lhs->targets.data(), (void **)targets));

      expansion_t expand{static_cast<expansion_t::contents_t *>(input->data),
                         input->bytes, input->n_digits};
      expansion->M_to_T(targets, &targets[lhs->targets.n()], input->scale);
      expansion->release();

      hpx_gas_unpin(lhs->targets.data());

      lhs->arrived += 1;
    } else if (*code == kLtoT) {
      LtoT *input = static_cast<LtoT *>(rhs);

      // The LCO data contains the reference to the targets, which must be
      // pinned.
      target_t *targets{nullptr}
      // NOTE: This should succeed because we place the LCO at the same locality
      // as the actual target data.
      assert(hpx_gas_try_pin(lhs->targets.data(), (void **)targets));

      expansion_t expand{static_cast<expansion_t::contents_t *>(input->data),
                         input->bytes, input->n_digits};
      expansion->L_to_T(lhs->targets, &lhs->targets[lhs->count], input->scale);
      expansion->release();

      hpx_gas_unpin(lhs->targets.data());

      lhs->arrived += 1;
    } else if (*code == kFinish) {
      lhs->finished = 1;
    } else {
      assert(0 && "Incorrect code to TargetLCO");
    }
  }

  static bool predicate_handler(TargetRefLCOData *i, size_t bytes) {
    return i->finished && (i->arrived == i->scheduled);
  }

  // The actual data for the object
  hpx_addr_t lco_;

  // the function identifiers
  static hpx_action_t init_;
  static hpx_action_t operation_;
  static hpx_action_t predicate_;
};


template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename, typename> class M>
hpx_action_t TargetLCO<S, T, E, M>::init_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename, typename> class M>
hpx_action_t TargetLCO<S, T, E, M>::operation_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename, typename> class M>
hpx_action_t TargetLCO<S, T, E, M>::predicate_ = HPX_ACTION_NULL;


} // namespace dashmm


#endif // __DASHMM_PARTICLE_H__
