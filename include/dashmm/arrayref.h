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


#ifndef __DASHMM_ARRAY_REF_H__
#define __DASHMM_ARRAY_REF_H__


/// \file
/// \brief ArrayRef object definition


#include <cstdlib>

#include <hpx/hpx.h>


namespace dashmm {


/// Reference to a set of records in a dashmm::Array object
///
/// This is a reference object, meaning that it refers to the Record data
/// but does not contain those data. As such, one can pass this
/// class by value without worry. Also, because this is a reference, this
/// object should not be used to destroy the memory to which it refers.
///
/// NOTE: This object is only intended to be used from an HPX-5 thread. Only
/// slice() actually bears this restriction, but the remaining routines are of
/// little use without access to slice().
///
/// This class is a template parameterized by the Record type.
template <typename Record>
class ArrayRef {
 public:
  using record_t = Record;
  using arrayref_t = ArrayRef<Record>;

  /// Default constructor.
  ArrayRef() : data_{nullptr}, n_{0} { }

  /// Construct from a specific address and counts.
  ArrayRef(Record *data, size_t n) : data_{data}, n_{n} { }

  /// Returns if the reference is valid
  bool valid() const {return data_ != nullptr;}

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
      return arrayref_t{nullptr, 0};
    }
    if (offset + n > n_) {
      return arrayref_t{nullptr, 0};
    }
    Record *addr = data_ + offset;
    return arrayref_t{addr, n};
  }

  /// Returns the number of Source records referred to.
  size_t n() const {return n_;}

  /// Returns the global address of the referred to data.
  Record *data() const {return data_;}

 private:
  /// Address of the first Record referred to by this object
  Record *data_;
  /// Number of Records referred to by this object
  size_t n_;
};


} // namespace dashmm


#endif // __DASHMM_SOURCE_REF_H__
