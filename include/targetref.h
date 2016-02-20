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


#ifndef __DASHMM_TARGET_REF_H__
#define __DASHMM_TARGET_REF_H__


/// \file include/targetref.h
/// \brief TargetRef object definition


#include <cstdlib>

#include <hpx/hpx.h>


namespace dashmm {


/// Target concept in DASHMM
///
/// To work as a Target type in DASHMM, a type must be trivially copyable and
/// must provide two (or more) members. The first is the location of the target
/// given by:
///
/// Point position
///
/// Second, one must provide a member for the output of the given expansion.
/// Often this is just
///
/// std::complex<double> phi
///
/// but it could be different for a given Expansion. See the documentation for
/// the Expansion of interest to determine the required members of Target for
/// using that Expansion. Some Expansion may produce multiple pieces of data
/// and so a target must supply each member.
///
/// Failure to provide the required members will result in compilation errors.


/// Reference to a set of targets
///
/// This is a reference object, meaning that it refers to the Target data in
/// the GAS, but does not contain those data. As such, one can pass this
/// class by value without worry. Also, because this is a reference, this object
/// cannot be used to destroy the memory in the global address space to which
/// it refers.
///
/// This class is a template parametrized by the Target type.
template <typename Target>
class TargetRef {
 public:
  using target_t = Target;
  using targetref_t = TargetRef<Target>;

  /// Default constructor.
  TargetRef()
      : data_{HPX_NULL}, n_{0}, n_tot_{0} { }

  /// Construct from a specific address and count.
  TargetRef(hpx_addr_t data, size_t n, size_t n_tot)
      : data_{data}, n_{n}, n_tot_{n_tot} { }

  /// Returns if the reference is valid
  bool valid() const {return data_ != HPX_NULL;}

  /// Get a reference to a slice of the current reference
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
  targetref_t slice(size_t offset, size_t n) const {
    if (offset > n_) {
      return targetref_t{HPX_NULL, 0, 0};
    }
    if (offset + n > n_) {
      return targetref_t{HPX_NULL, 0, 0};
    }
    hpx_addr_t addr = hpx_addr_add(data_, sizeof(Target) * offset,
                                   sizeof(Target) * n_tot_);
    return targetref_t{addr, n, n_tot_};
  }

  /// Returns the number of Source records referred to.
  size_t n() const {return n_;}

  /// Returns the total number of Source records in the underlying allocation.
  size_t n_tot() const {return n_tot_;}

  /// Returns the global address of the referred to data.
  hpx_addr_t data() const {return data_;}

 private:
  /// Address of the first Target referred to by this object
  hpx_addr_t data_;
  /// Number of Targets referred to by this object
  size_t n_;
  /// Total number of Targets in the underlying GAS allocation
  size_t n_tot_;
};


} // namespace dashmm


#endif // __DASHMM_TARGET_REF_H__
