// =============================================================================
//  This file is part of:
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  DASHMM is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  DASHMM is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with DASHMM. If not, see <http://www.gnu.org/licenses/>.
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

.
/// This class represents a pinned copy of a reference to an Array
///
/// An Array object represents data that lives in the global address space.
/// An ArrayRef object represents a portion of that array. An ArrayData
/// represents a locally accessible portion of the Array data.
///
/// This object is similar in some ways to a smart pointer, but instead of
/// managing the lifetime of the allocation, this will manage the lifetime of
/// an HPX-5 pin operation. To access the underlying data, one uses the
/// value() member.
template <typename Record>
class ArrayData {
 public:
  /// Default constructor
  ArrayData() : data_{HPX_NULL}, local_{nullptr} { }

  /// Construct from a global address
  explicit ArrayData(hpx_addr_t data) : data_{data}, local_{nullptr} {
    if (data_ != HPX_NULL) {
      assert(hpx_gas_try_pin(data_, (void **)&local_));
    }
  }

  /// Copy construction
  ArrayData(const ArrayData<Record> &other) {
    data_ = other.data_;
    if (data_ != HPX_NULL) {
      assert(hpx_gas_try_pin(data_, (void **)&local_));
    }
  }

  /// Move construction
  ArrayData(ArrayData<Record> &&other) {
    data_ = other.data_;
    local_ = other.local_;
    other.local_ = nullptr;
    other.data_ = HPX_NULL;
  }

  /// The destructor unpins the memory
  ~ArrayData() {
    if (data_ != HPX_NULL) {
      hpx_gas_unpin(data_);
    }
  }

  /// Copy assignment
  ArrayData<Record> &operator=(const ArrayData<Record> &other) {
    if (data_ != HPX_NULL) {
      hpx_gas_unpin(data_);
      local_ = nullptr;
    }

    data_ = other.data_;
    assert(hpx_gas_try_pin(data_, (void **)&local_));

    return *this;
  }

  /// Move assignment
  ArrayData<Record> &operator=(ArrayData<Record> &&other) {
    if (data_ != HPX_NULL) {
      hpx_gas_unpin(data_);
      local_ = nullptr;
    }

    data_ = other.data_;
    local_ = other.local_;

    other.data_ = HPX_NULL;
    other.local_ = nullptr;

    return *this;
  }

  /// Return the local address
  Record *value() const {return local_;}

 private:
  hpx_addr_t data_;
  Record *local_;
};


/// Reference to a set of records in a dashmm::Array object
///
/// This is a reference object, meaning that it refers to the Record data in
/// the GAS, but does not contain those data. As such, one can pass this
/// class by value without worry. Also, because this is a reference, this
/// object cannot be used to destroy the memory in the global address space to
/// which it refers.
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
  ArrayRef() : data_{HPX_NULL}, n_{0}, n_tot_{0} { }

  /// Construct from a specific address and counts.
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

  /// Return a local reference
  ArrayData<Record> pin() const {
    return ArrayData<Record>(data_);
  }

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
