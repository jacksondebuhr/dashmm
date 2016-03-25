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


#ifndef __DASHMM_ARRAY_REF_H__
#define __DASHMM_ARRAY_REF_H__


/// \file include/dashmm/arrayref.h
/// \brief ArrayRef object definition


#include <cstdlib>

#include <hpx/hpx.h>


namespace dashmm {


/// Reference to a set of records in a dashmm::Array object
///
/// This is a reference object, meaning that it refers to the Record data in
/// the GAS, but does not contain those data. As such, one can pass this
/// class by value without worry. Also, because this is a reference, this
/// object cannot be used to destroy the memory in the global address space to
/// which it refers.
///
/// This class is a template parameterized by the Record type.
template <typename Record>
class ArrayRef {
 public:
  using record_t = Record;
  using arrayref_t = ArrayRef<Record>;

  /// Default constructor.
  ArrayRef() : data_{HPX_NULL}, n_{0}, n_tot_{0} { }

  /// Construct from a specific address and count.
  ArrayRef(hpx_addr_t data, size_t n, size_t n_tot)
      : data_{data}, n_{n}, n_tot_{n_tot} { }

  /// Returns if the reference is valid
  bool valid() const {return data_ != HPX_NULL;}

  /// Get a reference to a slice of the current reference
  ///
  /// This will return an ArrayRef to a consecutive chunk of the records that
  /// begin at an offset from the start of this reference and that will contain
  /// n entries. If the input arguments to this method are invalid, then an
  /// invalid reference will be returned.
  ///
  /// \param offset - the offset from the start of this reference
  /// \param n - the number of entries in the slice
  ///
  /// \returns - the resulting ArrayRef; may be invalid.
  arrayref_t slice(size_t offset, size_t n) const {
    if (offset > n_) {
      return arrayref_t{HPX_NULL, 0, 0};
    }
    if (offset + n > n_) {
      return arrayref_t{HPX_NULL, 0, 0};
    }
    hpx_addr_t addr = hpx_addr_add(data_, sizeof(Record) * offset,
                                   sizeof(Record) * n_tot_);
    return arrayref_t{addr, n, n_tot_};
  }

  /// Returns the number of Source records referred to.
  size_t n() const {return n_;}

  /// Returns the total number of Source records in the underlying allocation.
  size_t n_tot() const {return n_tot_;}

  /// Returns the global address of the referred to data.
  hpx_addr_t data() const {return data_;}

 private:
  /// Address of the first Record referred to by this object
  hpx_addr_t data_;
  /// Number of Records referred to by this object
  size_t n_;
  /// Total number of Records in the underlying GAS allocation
  size_t n_tot_;
};


} // namespace dashmm


#endif // __DASHMM_SOURCE_REF_H__
