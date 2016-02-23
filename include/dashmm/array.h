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


#ifndef __DASHMM_ARRAY_H__
#define __DASHMM_ARRAY_H__


/// \file include/array.h
/// \brief Definitions needed to interact with DASHMM array objects.


#include <cstdlib>

#include <hpx/hpx.h>

#include "dashmm/types.h"


namespace dashmm {


/// Meta data for the array object
///
/// Array objects in DASHMM are two part objects in the global address space.
/// The ObjectHandle points to the meta data which is the information below.
/// Somewhere else in the global address space will be the data itself.
/// Currently, as DASHMM is specialized to SMP operation, the array data will
/// be on the same locality. In the future, this may change, and so we separate
/// the array meta data from the array itself, as the required metadata may
/// need to change as the data layout becomes more flexible.
struct ArrayMetaData {
  hpx_addr_t data;    /// global address of the array data
  size_t count;       /// the number of records in the array
  size_t size;        /// the size (bytes) of each record
};


/// Action for Array allocation
extern hpx_action_t allocate_array_action;

/// Action for Array deallocation
extern hpx_action_t deallocate_array_action;

/// Action for putting data into an Array
extern hpx_action_t array_put_action;

/// Action for getting data from an Array
extern hpx_action_t array_get_action;


/// Array class
///
/// Arrays are template classes parameterized over the data type of the
/// records stored in the array. The only requirement on the record type is
/// that it be trivally copyable.
///
/// This object is a reference to data in the global address space. As such,
/// is stores very little information, and so passing this object by value
/// does not have to be avoided.
///
/// Currently, this class should not be used inside an HPX-5 thread. Instead,
/// this is intended to be used from user code that makes calls to DASHMM.
template <typename T>
class Array {
 public:
  /// This creates the Array from the global address of an existing array
  ///
  /// Note that this address is the address of the Array meta data, and not
  /// the address of the records.
  ///
  /// \param data - the global address of the array meta data.
  Array(hpx_addr_t data = HPX_NULL) : data_{data} { }

  /// Returns if the Array is valid.
  ///
  /// An Array is valid if it refers to an array in the global address space.
  bool valid() const {return data != HPX_NULL;}

  /// Return the global address of the Array meta data.
  hpx_addr_t data() const {return data_;}

  /// This creates an array object
  ///
  /// \param record_count - the number of records that will be in the array.
  ReturnCode allocate(size_t record_count) {
    hpx_addr_t *dataout = &data_; //YIKES!
    int runcode;
    int *arg = &runcode;
    size_t size = sizeof(T);
    int err = hpx_run(&allocate_array_action, &record_count,
                      &size, &dataout, &arg);
    if (HPX_SUCCESS == err) {
      return kSuccess;
    } else {
      return static_cast<ReturnCode>(runcode);
    }
  }

  /// Destroy the Array
  ///
  /// This will destroy all global allocations associated with the Array.
  ///
  /// \returns - kSuccess on success; kRuntimeError otherwise.
  ReturnCode destroy() {
    if (HPX_SUCCESS == hpx_run(&deallocate_array_action, &data_)) {
      return kSuccess;
    } else {
      return kRuntimeError;
    }
  };

  /// Get data from an Array
  ///
  /// This retrieves data from the global address space and places it in
  /// 'normal' memory. The range of records to retrieve is specified using
  /// indices into the array, with the final index indicating one past the
  /// end of the range of interest.
  ///
  /// \param first - the first index of the range of interest
  /// \param last - one past the last index of the range of interest
  /// \param out_data - buffer into which the data will be placed
  ///
  /// \returns - kSuccess if successful; kRuntimeError if there is an error in
  ///            the runtime; kDomainError if the provided index range is
  ///            inconsistent with the Array object.
  ReturnCode get(size_t first, size_t last, T *out_data) {
    int runcode;
    int *arg = &runcode;
    hpx_run(&array_get_action, &data_, &first, &last, &arg, &out_data);
    return static_cast<ReturnCode>(runcode);
  }

  /// Put data into an Array
  ///
  /// This puts data from 'normal' memory into the specified range of records
  /// in the global address space. The range is specified with indices into the
  /// array, with the final index indicatign one past the end of the range of
  /// interest.
  ///
  /// \param first - the first index of the range of interest
  /// \param last - one past the last index of the range of interest
  /// \param in_data - buffer from which data will be copied
  ///
  /// \returns - kSuccess on success; kRuntimeError is there is an error in the
  ///            runtime; kDomainError if the provided index range is
  ///            inconsistent with the Array object.
  ReturnCode put(size_t first, size_t last, T *in_data) {
    int runcode;
    int *arg = &runcode;
    hpx_run(&array_put_action, &data_, &first, &last, &arg, &in_data);
    return static_cast<ReturnCode>(runcode);
  }

 private:
  /// The global address of the Array meta data.
  hpx_addr_t data_;
};



} // namespace dashmm


#endif // __DASHMM_ARRAY_H__
