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


/// \file include/dashmm/array.h
/// \brief Definitions needed to interact with DASHMM array objects.


#include <cstdlib>

#include <hpx/hpx.h>

#include "dashmm/arraymetadata.h"
#include "dashmm/arraymapaction.h"
#include "dashmm/arrayref.h"
#include "dashmm/types.h"


namespace dashmm {


/// Action for Array meta data allocation
extern hpx_action_t allocate_array_meta_action;

/// Action for local portions of Array allocation
extern hpx_action_t allocate_local_work_action;

/// Action to delete organizational stucture for Array allocation
extern hpx_action_t allocate_array_destroy_reducer_action;

/// Action for Array deallocation
extern hpx_action_t deallocate_array_action;

/// Action for putting data into an Array
extern hpx_action_t array_put_action;

/// Action for getting data from an Array
extern hpx_action_t array_get_action;

/// Action for getting local counts
extern hpx_action_t array_local_count_action;

/// Action for getting total count
extern hpx_action_t array_total_count_action;


struct ArrayMetaAllocRunReturn {
  hpx_addr_t meta;
  hpx_addr_t reducer;
  int code;
};


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
  /// NOTE: This only tests that this object points to someplace in the global
  /// address space. This is an incomplete test. A more complete test would
  /// examine the data stored at the specified address for compatibility.
  ///
  /// \returns - true if this object refers to non-null global memory.
  bool valid() const {return data != HPX_NULL;}

  /// Return the global address of the Array meta data.
  ///
  /// \returns - the global address of the Array meta data.
  hpx_addr_t data() const {return data_;}

  /// Return the number of records in the local portion of this object
  ///
  /// This is a SPMD style operation; it is collective, but each rank will
  /// receive a different value giving the local portion on the calling rank.
  ///
  /// \returns - the number of records owned by this rank,
  size_t count() const {
    size_t retval{0};
    hpx_run_spmd(&array_local_count_action, &retval, &data_);
    return retval;
  }

  /// Return the number of records in the entire array
  ///
  /// This is a SPMD style operation; it is collective, and all ranks will
  /// receive the same result.
  ///
  /// \returns - the total number of records in all ranks.
  size_t length() const {
    size_t retval{0};
    hpx_run_spmd(&array_total_count_action, &retval, &data_);
    return retval;
  }

  /// This creates an array object
  ///
  /// This is called from the SPMD user-application. This cannot be used
  /// inside an HPX-5 thread. Each rank of the user application will provide
  /// its own record count. The resulting array will be created with the
  /// provided distribution of records. One or more of these counts can be
  /// zero provided there is at least one non-zero count from some rank.
  ///
  /// \param record_count - the number of records that will be in the array
  ///                       for this locality.
  ///
  /// \returns - kDomainError if the object already has an allocation;
  ///            kAllocationError if the global memory cannot be allocated;
  ///            kRuntimeError if there is an error in the runtime; or
  ///            kSuccess otherwise.
  ReturnCode allocate(size_t record_count) {
    if (data_ != HPX_NULL) {
      // If the object already has data, do not allocate new data.
      return kDomainError;
    }

    ArrayMetaAllocRunReturn metadata{HPX_NULL, HPX_NULL, 0};
    hpx_run(&allocate_array_meta_action, &metadata, nullptr, 0);
    if (metadata.code != kSuccess) {
      return static_cast<ReturnCode>(metadata.code);
    }
    data_ = metadata.meta;

    ReturnCode retval{kSuccess};

    size_t size = sizeof(T);
    int runcode{0};
    hpx_run_spmd(&allocate_local_work_action, &runcode,
                 &metadata.meta, &metadata.reducer, &size, &record_count);
    if (runcode != kSuccess) {
      // NOTE: If this branch is triggered, only some ranks will do this.
      //  so there is an issue to sort out about cleaning up the resulting
      //  mess.
      // We can likely just get another reduction, and send it back, and
      //  reduce on it before we leave the local work action. Then everyone
      //  would agree on the result.
      // TODO; probably do that...
      retval = kAllocationError;
    }

    hpx_run(&allocate_array_destroy_reducer_action, nullptr, &metadata.reducer);

    return retval;
  }

  /// Destroy the Array
  ///
  /// This will destroy all global allocations associated with the Array.
  /// If this does not return kSuccess, the global memory will not have been
  /// released.
  ///
  /// \returns - kSuccess on success; kRuntimeError otherwise.
  ReturnCode destroy() {
    if (HPX_SUCCESS == hpx_run(&deallocate_array_action, nullptr, &data_)) {
      data_ = HPX_NULL;
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
  /// This routine is performed on a rank-by-rank basis, meaning that the
  /// input arguments can be different from rank to rank. Also, this routine
  /// will only get data from the portion of the array that is local to the
  /// rank. The indices given are relative to the start of the local portion.
  ///
  /// \param first - the first index of the range of interest
  /// \param last - one past the last index of the range of interest
  /// \param out_data - buffer into which the data will be placed
  ///
  /// \returns - kSuccess if successful; kRuntimeError if there is an error in
  ///            the runtime; kDomainError if the provided index range is
  ///            inconsistent with the Array object.
  ReturnCode get(size_t first, size_t last, T *out_data) {
    int runcode = kSuccess;
    // TODO: the error handling with this is now strange. Each rank can
    // have a different runcode here...
    hpx_run_spmd(&array_get_action, &runcode, &data_, &first, &last, &out_data);
    return static_cast<ReturnCode>(runcode);
  }

  /// Put data into an Array
  ///
  /// This puts data from 'normal' memory into the specified range of records
  /// in the global address space. The range is specified with indices into the
  /// array, with the final index indicating one past the end of the range of
  /// interest.
  ///
  /// This routine is performed on a rank-by-rank basis, meaning that the
  /// input arguments can be different from rank to rank. Also, this routine
  /// will only put data into the portion of the array that is local to the
  /// rank. The indices given are relative to the start of the local portion.
  ///
  /// \param first - the first index of the range of interest
  /// \param last - one past the last index of the range of interest
  /// \param in_data - buffer from which data will be copied
  ///
  /// \returns - kSuccess on success; kRuntimeError is there is an error in the
  ///            runtime; kDomainError if the provided index range is
  ///            inconsistent with the Array object.
  ReturnCode put(size_t first, size_t last, T *in_data) {
    int runcode = kSuccess;
    // TODO: the error handling with this is now strange. Each rank can
    // have a different runcode here...
    hpx_run_spmd(&array_put_action, &runcode, &data_, &first, &last, &in_data);
    return static_cast<ReturnCode>(runcode);
  }

  /// Produce an ArrayRef from the Array
  ///
  /// Largely speaking, this will not be needed by DASHMM users. This array
  /// reference will refer to the entirety of the local part at the given
  /// locality
  ///
  /// NOTE: This routine cannot be called outside of an HPX-5 thread.
  ///
  /// \param rank - the rank to get a reference to. Provide -1 to indicate
  ///               whatever is the caller's rank.
  ///
  /// \returns - an ArrayRef indicating the entirety of this array.
  ArrayRef<T> ref(int rank = -1) const {
    if (rank == -1) {
      rank = hpx_get_my_rank();
    }
    hpx_addr_t global = hpx_addr_add(data_,
                                     sizeof(ArrayMetaData) * rank,
                                     sizeof(ArrayMetaData));
    ArrayMetaData *local{nullptr};
    assert(hpx_gas_try_pin(global, (void **)&local));
    ArrayRef<T> retval{local->data, local->local_count, local->local_count};
    hpx_gas_unpin(global);
    return retval;
  }

  /// Map an action onto each record in the Array
  ///
  /// This will cause the action represented by the given argument, to be
  /// invoked on all the entries of the array. The action will ultimatly work
  /// on segments of the array. The environment is provided unmodified to each
  /// segment. Please see the ArrayMapAction for more details.
  ///
  /// This acts on a rank-by-rank basis. Each rank will participate, and will
  /// handle the records in the array owned by the rank. In principle, the
  /// environment passed to each rank might have a different value.
  ///
  /// \param act - the action to perform on the entries of the array.
  /// \param env - the environment to use in the action.
  ///
  /// \return - kSuccess
  template <typename E>
  ReturnCode map(const ArrayMapAction<T, E> &act, const E *env) {
    hpx_run_spmd(&act.root_, nullptr, &act.leaf_, &env, &data_);
    return kSuccess;
  }

 private:
  /// The global address of the Array meta data.
  hpx_addr_t data_;
};


} // namespace dashmm


#endif // __DASHMM_ARRAY_H__
