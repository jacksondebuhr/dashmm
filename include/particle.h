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
  /// this oject only destroys the reference.
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


/// The data needed for target locations.
struct Target {
  Point position;
  dcomplex_t phi;
  /// DASHMM rearranges the input locations so we store the index from the
  /// original data so that we can copy results into the correct place.
  size_t index;
};


/// Reference to target data
///
/// This object is a reference, but unlike SourceRef, it refers to an LCO.
/// The LCO is a user-defined LCO that handles the concurrent contribution to
/// the potential of the target points. The data of the LCO is the actual
/// Target data and a few additional pieces of management data.
class TargetRef {
 public:
  /// Construct a default object
  TargetRef() : data_{HPX_NULL}, n_{0} { }

  /// Construct a reference from an existing global address
  TargetRef(hpx_addr_t data, int n) : data_{data}, n_{n} { }

  /// Construct a reference from input targets. This will create the LCO
  TargetRef(Target *targets, int n);

  /// As a reference, destruction of the referred to object must be done
  /// via destroy.
  void destroy();

  /// The number of target points represented by the object
  int n() const {return n_;}

  /// The global address of the referred to object
  hpx_addr_t data() const {return data_;}

  /// Indicate to the underlying LCO that all operations have been scheduled
  void finalize() const;

  /// Indicate to the underlying LCO that it should expect \param num more
  /// operations
  void schedule(int num) const;

  /// Contribute a S->T operation to the referred targets
  ///
  /// \param type - the type of the expansion serving the S->T operation
  /// \param n - the number of sources
  /// \param sources - the sources themselves
  void contribute_S_to_T(int type, int n, Source *sources) const;

  /// Contribute a M->T operation to the referred targets
  ///
  /// \param type - the type of expansion serving the M->T operation
  /// \param bytes - the size of the serialized expansion data
  /// \param data - the serialized expansion data
  /// \param n_digits - accuracy of the expansion
  /// \param scale - scaling factor
  void contribute_M_to_T(int type, size_t bytes, void *data,
                         int n_digits, double scale) const;

  /// Contribute a L->T operation to the referred targets
  ///
  /// \param type - the type of expansion serving the L->T operation
  /// \param bytes - the size of the serialized expansion data
  /// \param data - the serialized expansion data
  /// \param n_digits - accuracy of the expansion
  /// \param scale - scaling factor
  void contribute_L_to_T(int type, size_t bytes, void *data,
                         int n_digits, double scale) const;

 private:
  hpx_addr_t data_;
  int n_;
};


} // namespace dashmm


#endif // __DASHMM_PARTICLE_H__
