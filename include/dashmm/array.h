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


#ifndef __DASHMM_ARRAY_H__
#define __DASHMM_ARRAY_H__


/// \file
/// \brief Definitions needed to interact with DASHMM array objects.

#include <cassert>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <memory>

#include <hpx/hpx.h>

#include "dashmm/arraymetadata.h"
#include "dashmm/arraymapaction.h"
#include "dashmm/arrayref.h"
#include "dashmm/reductionops.h"
#include "dashmm/types.h"

#include "builtins/trivialserializer.h"


namespace dashmm {


template <typename T>
class ArrayRegistrar;


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
  bool valid() const {return data_ != HPX_NULL;}

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
    assert(valid());
    size_t retval{0};
    hpx_run_spmd(&array_local_count_, &retval, &data_);
    return retval;
  }

  /// Return the number of records in the entire array
  ///
  /// This is a SPMD style operation; it is collective, and all ranks will
  /// receive the same result.
  ///
  /// \returns - the total number of records in all ranks.
  size_t length() const {
    assert(valid());
    size_t retval{0};
    hpx_run_spmd(&array_total_count_, &retval, &data_);
    return retval;
  }

  /// This allocates an array object
  ///
  /// This is called from the SPMD user-application. This cannot be used
  /// inside an HPX-5 thread. Each rank of the user application will provide
  /// its own @p record count. The resulting array will be created with the
  /// provided distribution of records. One or more of these counts can be
  /// zero provided there is at least one non-zero count from some rank.
  ///
  /// Optionally, memory that has already been allocated can be provided as
  /// a second argument to this function, which will then be used as the
  /// local segment of the array for the calling rank. DASHMM will assume
  /// ownership of a segment provided in this fashion.  Each rank can provide
  /// a segment and the resulting Array will be formed of the provided
  /// segments. If this options is omitted or is nullptr, DASHMM will allocate
  /// memory, which must be accessed via other methods.
  ///
  /// \param record_count - the number of records that will be in the array
  ///                       for this locality.
  /// \param segment - some previously allocated memory to use as the segment
  ///                  of the array for this rank.
  ///
  /// \returns - kDomainError if the object already has an allocation;
  ///            kAllocationError if the global memory cannot be allocated;
  ///            kRuntimeError if there is an error in the runtime; or
  ///            kSuccess otherwise.
  ReturnCode allocate(size_t record_count, T *segment = nullptr) {
    if (data_ != HPX_NULL) {
      // If the object already has data, do not allocate new data.
      return kDomainError;
    }

    ArrayMetaAllocRunReturn metadata{HPX_NULL, HPX_NULL, HPX_NULL, 0};
    hpx_run(&allocate_array_meta_, &metadata);
    if (metadata.code != kSuccess) {
      return static_cast<ReturnCode>(metadata.code);
    }
    data_ = metadata.meta;

    ReturnCode retval{kSuccess};

    int runcode{0};
    hpx_run_spmd(&allocate_local_work_, &runcode,
                 &metadata.meta, &metadata.reducer, &metadata.retcode,
                 &record_count, &segment);
    if (runcode != kSuccess) {
      retval = kAllocationError;
    }

    hpx_run(&allocate_array_destroy_reducer_, nullptr,
            &metadata.reducer, &metadata.retcode);

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
    if (HPX_SUCCESS == hpx_run(&deallocate_array_, nullptr, &data_)) {
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
    assert(valid());

    hpx_addr_t reducer{HPX_NULL};
    hpx_run(&get_or_put_retcode_reducer_, &reducer);

    int runcode{kSuccess};
    hpx_run_spmd(&array_get_, &runcode, &data_, &first, &last, &out_data,
                 &reducer);

    hpx_run(&get_or_put_reducer_delete_, nullptr, &reducer);

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
    assert(valid());

    hpx_addr_t reducer{HPX_NULL};
    hpx_run(&get_or_put_retcode_reducer_, &reducer);

    int runcode = kSuccess;
    hpx_run_spmd(&array_put_, &runcode, &data_, &first, &last, &in_data,
                 &reducer);

    hpx_run(&get_or_put_reducer_delete_, nullptr, &reducer);

    return static_cast<ReturnCode>(runcode);
  }

  /// Get access to the local segment
  ///
  /// This will return the value of segment to be the address of the local
  /// segment of the array on this rank. This is a collective call, and so
  /// each rank will receive the address of its segment. Further, as it is
  /// possible that other DASHMM routines will sort and rearrange Array data,
  /// the length of the segment is also returned.
  ///
  /// Calls to other DASHMM routines may invalidate the address returned by
  /// this routine.
  ///
  /// \param count [out] - length of segment
  ///
  /// \returns - address of local segment; this can be nullptr when the segment
  ///            is empty.
  T *segment(size_t &count) {
    SegmentReturn retval{nullptr, 0};
    hpx_run_spmd(&array_segment_request_, &retval, &data_);

    count = retval.count;
    return retval.segment;
  }

  /// Produce an ArrayRef from the Array
  ///
  /// Largely speaking, this will not be needed by DASHMM users. This array
  /// reference will refer to the entirety of the local part at the given
  /// locality
  ///
  /// NOTE: This routine cannot be called outside of an HPX-5 thread.
  ///
  /// \returns - an ArrayRef indicating the entirety of this array.
  ArrayRef<T> ref() const {
    int rank = hpx_get_my_rank();
    hpx_addr_t global = hpx_addr_add(data_,
                                     sizeof(ArrayMetaData<T>) * rank,
                                     sizeof(ArrayMetaData<T>));
    ArrayMetaData<T> *local{nullptr};
    assert(hpx_gas_try_pin(global, (void **)&local));
    ArrayRef<T> retval{local->data, local->local_count};
    hpx_gas_unpin(global);
    return retval;
  }

  /// Replace the local segment. If these change the overall total number of
  /// records, that will have to be fixed with a call to resum().
  ///
  /// NOTE: This can only be called from inside an HPX-5 thread, and as such,
  /// It is suggested that the casual user not use this routine.
  ///
  /// \param ref - an ArrayRef giving the new segment to install for this rank
  ///
  /// \returns - the address of the previous data segment; the user assumes
  ///            ownership of this data.
  T *replace(ArrayRef<T> ref) {
    int rank = hpx_get_my_rank();
    hpx_addr_t global = hpx_addr_add(data_,
                                     sizeof(ArrayMetaData<T>) * rank,
                                     sizeof(ArrayMetaData<T>));
    ArrayMetaData<T> *local{nullptr};
    assert(hpx_gas_try_pin(global, (void **)&local));

    T *retval = local->data;
    local->data = ref.data();
    local->local_count = ref.n();

    hpx_gas_unpin(global);
    return retval;
  }

  /// Map an action onto each record in the Array
  ///
  /// This will cause the action represented by @p act, to be
  /// invoked on all the entries of the array. The action will ultimatly work
  /// on segments of the array. The environment. @p env, is provided unmodified
  /// to each segment. Please see the ArrayMapAction for more details.
  ///
  /// This acts on a rank-by-rank basis. Each rank will participate, and will
  /// handle the records in the array owned by the rank. In principle, the
  /// environment passed to each rank might have a different value.
  ///
  /// \param act - the action to perform on the entries of the array.
  /// \param env - the environment to use in the action.
  ///
  /// \return - kSuccess
  template <typename E, int F>
  ReturnCode map(const ArrayMapAction<T, E, F> &act, const E *env) {
    assert(valid());
    hpx_run_spmd(&act.root_, nullptr, &act.leaf_, &env, &data_);
    return kSuccess;
  }

  /// Collect all of an array's data into a single local array
  ///
  /// This will return a newly allocated array containing all of the records
  /// represented by this Array (note the difference in spelling). The
  /// ordering of the records is not guaranteed, and so users should track
  /// record identity in some way. This does not modify the data in the
  /// global address space. Instead, this copies the entire array into one
  /// local segment. This is predominantly intended as an ease-of-use method,
  /// and should not be expected to perform well.
  ///
  /// Only rank zero will return data. All other ranks will receive nullptr
  /// from this routine.
  ///
  /// \returns - smart pointer containing the local allocation
  std::unique_ptr<T[]> collect() {
    assert(valid());

    // I need a step to put the serializer into a known location...

    T *retval{nullptr};
    hpx_run(&array_collect_, &retval, &data_);

    if (hpx_get_my_rank()) {
      retval = nullptr;
    }

    return std::unique_ptr<T[]>(retval);
  }

  /// Set the serialization manager for the Array
  ///
  /// This will associate the given serialization manager object with this
  /// Array. This will then be used inside any DASHMM routine that needs it
  /// when dealing with this data.
  ///
  /// It is important that the same type be given to each rank when calling
  /// this method, otherwise inconsistencies could result.
  ///
  /// \param manager - the Serializer object to associate with this Array
  ///
  /// \returns - kSuccess on success; kRuntimeError if there was a problem in
  ///            the runtime.
  ReturnCode set_manager(std::unique_ptr<Serializer> manager) {
    Serializer *ptr = manager.release();
    if (HPX_SUCCESS != hpx_run_spmd(&array_set_manager_, nullptr,
                                    &ptr, &data_)) {
      return kRuntimeError;
    }
    return kSuccess;
  }

  /// Get serialization manager
  ///
  /// This returns the address of this rank's serialization manager of this
  /// Array.
  ///
  /// NOTE: This should only be called from inside an HPX-5 epoch, and should
  /// not be considered to be part of the user interface.
  ///
  /// \returns - the serializaion manager
  Serializer *get_manager() const {
    int rank = hpx_get_my_rank();
    hpx_addr_t global = hpx_addr_add(data_,
                                     sizeof(ArrayMetaData<T>) * rank,
                                     sizeof(ArrayMetaData<T>));
    ArrayMetaData<T> *local{nullptr};
    assert(hpx_gas_try_pin(global, (void **)&local));

    auto retval = local->manager;

    hpx_gas_unpin(global);
    return retval;
  }

 private:
  friend class ArrayRegistrar<T>;

  hpx_addr_t data_;

  static hpx_action_t array_local_count_;
  static hpx_action_t array_total_count_;
  static hpx_action_t allocate_array_meta_;
  static hpx_action_t allocate_local_work_;
  static hpx_action_t allocate_array_destroy_reducer_;
  static hpx_action_t deallocate_array_;
  static hpx_action_t deallocate_array_local_;
  static hpx_action_t get_or_put_retcode_reducer_;
  static hpx_action_t get_or_put_reducer_delete_;
  static hpx_action_t array_get_;
  static hpx_action_t array_put_;
  static hpx_action_t array_segment_request_;
  static hpx_action_t array_collect_;
  static hpx_action_t array_collect_request_;
  static hpx_action_t array_collect_receive_;
  static hpx_action_t array_set_manager_;

  /// Return data from allocation action
  struct ArrayMetaAllocRunReturn {
    hpx_addr_t meta;
    hpx_addr_t reducer;
    hpx_addr_t retcode;
    int code;
  };

  /// Return data from segment request
  struct SegmentReturn {
    T *segment;
    size_t count;
  };


  /// Action to return the local length of the Array
  ///
  /// \param data - global address of the Array's meta data
  ///
  /// \returns - HPX_SUCCESS
  static int array_local_count_handler(hpx_addr_t data) {
    hpx_addr_t global = hpx_addr_add(data,
        sizeof(ArrayMetaData<T>) * hpx_get_my_rank(),
        sizeof(ArrayMetaData<T>));
    ArrayMetaData<T> *local{nullptr};
    assert(hpx_gas_try_pin(global, (void **)&local));

    size_t retval = local->local_count;

    hpx_gas_unpin(global);
    hpx_exit(sizeof(retval), &retval);
  }

  /// Action to return the Array's total length
  ///
  /// \param data - global address of the Array's meta data
  ///
  /// \returns - HPX_SUCCESS
  static int array_total_count_handler(hpx_addr_t data) {
    hpx_addr_t global = hpx_addr_add(data,
            sizeof(ArrayMetaData<T>) * hpx_get_my_rank(),
            sizeof(ArrayMetaData<T>));
    ArrayMetaData<T> *local{nullptr};
    assert(hpx_gas_try_pin(global, (void **)&local));

    size_t retval = local->total_count;

    hpx_gas_unpin(global);
    hpx_exit(sizeof(retval), &retval);
  }

  /// Action for allocating Array metadata
  ///
  /// This action is called on a single rank to allocate some collective
  /// data that is needed to allocate an array. The returned information is
  /// broadcast back to every rank during hpx_exit().
  ///
  /// \returns - HPX_SUCCESS
  static int allocate_array_meta_handler() {
    ArrayMetaAllocRunReturn retval{HPX_NULL, HPX_NULL, kSuccess};

    int ranks = hpx_get_num_ranks();
    retval.meta = hpx_gas_alloc_cyclic(ranks, sizeof(ArrayMetaData<T>), 0);
    if (retval.meta == HPX_NULL) {
      retval.code = kAllocationError;
      hpx_exit(sizeof(retval), &retval);
    }

    retval.reducer = hpx_lco_reduce_new(ranks, sizeof(size_t) * ranks,
                                        size_sum_ident, size_sum_op);
    if (retval.reducer == HPX_NULL) {
      hpx_gas_free_sync(retval.meta);
      retval.meta = HPX_NULL;
      retval.code = kAllocationError;
    }

    retval.retcode = hpx_lco_reduce_new(ranks, sizeof(int),
                                        int_max_ident_op, int_max_op);
    if (retval.retcode == HPX_NULL) {
      hpx_gas_free_sync(retval.meta);
      hpx_lco_delete_sync(retval.reducer);
      retval.meta = HPX_NULL;
      retval.reducer = HPX_NULL;
      retval.code = kAllocationError;
    }

    hpx_exit(sizeof(retval), &retval);
  }

  /// Action for allocation rank-local segments of an Array
  ///
  /// This action, called in a SPMD fashion on all ranks, will both participate
  /// in a reduction to count the total size of the array and will allocate the
  /// local segment when that reduction is complete.
  ///
  /// \param data - the global address of the ArrayMetaData
  /// \param reducer - the address of the reduction LCO used to compute the
  ///                   total array length
  /// \param record_size - the size of each record
  /// \param record_count - the number of records for this rank
  ///
  /// \returns - HPX_SUCCESS
  static int allocate_local_work_handler(hpx_addr_t data, hpx_addr_t reducer,
                                         hpx_addr_t retcode,
                                         size_t record_count,
                                         T *segment) {
    int ranks = hpx_get_num_ranks();
    int my_rank = hpx_get_my_rank();
    size_t *contrib = new size_t [ranks];
    assert(contrib);
    for (int i = 0; i < ranks; ++i) {
      if (i >= my_rank) {
        contrib[i] = record_count;
      } else {
        contrib[i] = 0;
      }
    }

    hpx_lco_set(reducer, sizeof(size_t) * ranks, contrib, HPX_NULL, HPX_NULL);
    hpx_lco_get(reducer, sizeof(size_t) * ranks, contrib);

    // Pin local part of the meta data
    hpx_addr_t global = hpx_addr_add(data,
        hpx_get_my_rank() * sizeof(ArrayMetaData<T>),
        sizeof(ArrayMetaData<T>));
    ArrayMetaData<T> *local{nullptr};
    assert(hpx_gas_try_pin(global, (void **)&local));
    local->local_count = record_count;
    local->total_count = contrib[ranks - 1];
    local->data = segment;
    local->manager = new TrivialSerializer<T>{};

    delete [] contrib;

    int retval{0};
    if (record_count && local->data == nullptr) {
      try {
        local->data = new T[record_count]{};
      } catch (std::bad_alloc &ba) {
        retval = kAllocationError;
      }
    }

    hpx_lco_set_lsync(retcode, sizeof(int), &retval, HPX_NULL);
    hpx_lco_get(retcode, sizeof(int), &retval);

    hpx_gas_unpin(global);
    hpx_exit(sizeof(retval), &retval);
  }

  /// Action to deallocate the reduction LCO used during Array allocation
  ///
  /// \param reducer - the global address of the reduction LCO
  ///
  /// \returns - HPX_SUCCESS
  static int allocate_array_destroy_reducer_handler(hpx_addr_t reducer,
                                                    hpx_addr_t retcode) {
    hpx_lco_delete_sync(reducer);
    hpx_lco_delete_sync(retcode);
    hpx_exit(0, nullptr);
  }

  /// Action that deallocates an array object
  ///
  /// This will delete the ArrayMetaData and the records themselves.
  ///
  /// \param obj - the global address of the array's meta data
  ///
  /// \returns - HPX_SUCCESS
  static int deallocate_array_handler(hpx_addr_t obj) {
    hpx_bcast_rsync(deallocate_array_local_, &obj);
    hpx_gas_free_sync(obj);
    hpx_exit(0, nullptr);
  }

  /// Action that deletes the local portion of an Array
  ///
  /// This action is the target of a broadcast. It will delete the local portion
  /// of the Array's global memory.
  ///
  /// \param meta - the global address of the Array's meta data.
  ///
  /// \returns - HPX_SUCCESS
  static int deallocate_array_local_handler(hpx_addr_t meta) {
    hpx_addr_t global = hpx_addr_add(meta,
          sizeof(ArrayMetaData<T>) * hpx_get_my_rank(),
          sizeof(ArrayMetaData<T>));
    ArrayMetaData<T> *local{nullptr};
    assert(hpx_gas_try_pin(global, (void **)&local));

    if (local->data != HPX_NULL) {
      delete [] local->data;
    }

    hpx_gas_unpin(global);

    return HPX_SUCCESS;
  }

  /// Action that allocates a reducer for return codes
  ///
  /// To be very careful, we collect the returns codes from each rank and
  /// take the maximum value.
  ///
  /// \returns - HPX_SUCCESS
  static int get_or_put_retcode_reducer_handler(void) {
    hpx_addr_t retval = hpx_lco_reduce_new(hpx_get_num_ranks(),
        sizeof(int), int_max_ident_op, int_max_op);
    hpx_exit(sizeof(retval), &retval);
  }

  /// Action that deletes a reducer for return codes
  ///
  /// \param lco - the LCO's address
  ///
  /// \returns - HPX_SUCCESS
  static int get_or_put_reducer_delete_handler(hpx_addr_t lco) {
    hpx_lco_delete_sync(lco);
    hpx_exit(0, nullptr);
  }

  /// Get data from an array object
  ///
  /// This will read records from a global array into @p out_data. This
  /// action is not well behaved as regards bad input arguments. This will
  /// likely be updated in the future.
  ///
  /// \param obj - the global address of the ArrayMetaData
  /// \param first - the first record to read
  /// \param last - the last (exclusive) record to read
  /// \param out_data - a buffer into which the read data is stored
  /// \param reducer - the LCO that will reduce the error codes from each rank
  ///
  /// \returns - HPX_SUCCESS
  static int array_get_handler(hpx_addr_t obj, size_t first, size_t last,
                               T *out_data, hpx_addr_t reducer) {
    int retval = kSuccess;

    hpx_addr_t global = hpx_addr_add(obj,
          sizeof(ArrayMetaData<T>) * hpx_get_my_rank(),
          sizeof(ArrayMetaData<T>));
    ArrayMetaData<T> *local{nullptr};
    if (!hpx_gas_try_pin(global, (void **)&local)) {
      retval = kRuntimeError;
    } else {
      if (last > local->local_count || last < first) {
        retval = kDomainError;
      } else if (local->local_count) {
        assert(local->data != nullptr);
        std::copy(&local->data[first], &local->data[last], out_data);
      }

      hpx_gas_unpin(global);
    }

    // reduce on error condition
    hpx_lco_set_lsync(reducer, sizeof(int), &retval, HPX_NULL);
    hpx_lco_get(reducer, sizeof(int), &retval);

    hpx_exit(sizeof(retval), &retval);
  }

  /// Put data into an array object
  ///
  /// This will fill records in a global array with @p in_data. This
  /// action is not well behaved as regards bad input arguments. This will likely
  /// be updated in the future.
  ///
  /// \param obj - the global address of the ArrayMetaData
  /// \param first - the first record to put data into
  /// \param last - the last (exclusive) record to put data into
  /// \param in_data - the data to copy into the records
  /// \param reducer - the LCO that will combine error codes from each rank
  ///
  /// \returns - HPX_SUCCESS
  static int array_put_handler(hpx_addr_t obj, size_t first, size_t last,
                               T *in_data, hpx_addr_t reducer) {
    int retval = kSuccess;

    hpx_addr_t global = hpx_addr_add(obj,
            sizeof(ArrayMetaData<T>) * hpx_get_my_rank(),
            sizeof(ArrayMetaData<T>));
    ArrayMetaData<T> *local{nullptr};
    if (!hpx_gas_try_pin(global, (void **)&local)) {
      retval = kRuntimeError;
    } else {
      if (last > local->local_count || last < first) {
        retval = kDomainError;
      } else if (local->local_count) {
        assert(local->data != nullptr);
        std::copy(in_data, &in_data[last], &local->data[first]);
      }

      hpx_gas_unpin(global);
    }

    hpx_lco_set_lsync(reducer, sizeof(int), &retval, HPX_NULL);
    hpx_lco_get(reducer, sizeof(int), &retval);

    hpx_exit(sizeof(retval), &retval);
  }

  static int array_segment_request_handler(hpx_addr_t data) {
    hpx_addr_t global = hpx_addr_add(data,
          sizeof(ArrayMetaData<T>) * hpx_get_my_rank(),
          sizeof(ArrayMetaData<T>));
    ArrayMetaData<T> *local{nullptr};
    assert(hpx_gas_try_pin(global, (void **)&local));

    // collect info
    SegmentReturn retval{};
    retval.segment = local->data;
    retval.count = local->local_count;

    hpx_gas_unpin(global);

    hpx_exit(sizeof(retval), &retval);
  }

  static int array_collect_handler(hpx_addr_t data) {
    // collect counts from GAS
    int my_rank = hpx_get_my_rank();
    int n_rank = hpx_get_num_ranks();

    ArrayMetaData<T> *meta = (ArrayMetaData<T> *)hpx_malloc_registered(
                                          sizeof(ArrayMetaData<T>) * n_rank);
    for (int i = 0; i < n_rank; ++i) {
      hpx_addr_t target = hpx_addr_add(data,
                                       sizeof(ArrayMetaData<T>) * i,
                                       sizeof(ArrayMetaData<T>));
      hpx_gas_memget_sync(&meta[i], target, sizeof(ArrayMetaData<T>));
    }

    // do some arithmetic
    size_t *offsets = new size_t[n_rank];
    offsets[0] = 0;
    for (int i = 1; i < n_rank; ++i) {
      offsets[i] = offsets[i - 1] + meta[i - 1].local_count;
    }

    // create AND gate for completion detection
    hpx_addr_t done = hpx_lco_and_new(n_rank);

    // allocate space for result
    T *retval{nullptr};
    retval = new T[meta[my_rank].total_count];

    // call on each rank with location and and gate and data
    for (int i = 0; i < n_rank; ++i) {
      if (i != my_rank) {
        T *location = retval + offsets[i];
        hpx_call(HPX_THERE(i), array_collect_request_, HPX_NULL,
                 &data, &done, &location);
      }
    }

    // copy my segment over
    {
      T *c_start = meta[my_rank].data;
      T *c_end = &meta[my_rank].data[meta[my_rank].local_count];
      T *d_start = &retval[offsets[my_rank]];
      std::copy(c_start, c_end, d_start);
    }
    hpx_lco_and_set(done, HPX_NULL);

    // wait for results
    hpx_lco_wait(done);
    hpx_lco_delete_sync(done);

    hpx_free_registered(meta);
    delete [] offsets;

    hpx_exit(sizeof(retval), &retval);
  }

  static int array_collect_request_handler(hpx_addr_t data, hpx_addr_t done,
                                           T *location) {
    hpx_addr_t global = hpx_addr_add(data,
              sizeof(ArrayMetaData<T>) * hpx_get_my_rank(),
              sizeof(ArrayMetaData<T>));
    ArrayMetaData<T> *local{nullptr};
    assert(hpx_gas_try_pin(global, (void **)&local));

    // get a parcel of the right size
    size_t arrsize{0};
    for (size_t i = 0; i < local->local_count; ++i) {
      // We have to count each individually because each object may have a
      // different serialized size.
      arrsize += local->manager->size(local->data + i);
    }
    size_t msgsize = sizeof(hpx_addr_t) + sizeof(T *) + arrsize;
    hpx_parcel_t *p = hpx_parcel_acquire(nullptr, msgsize);
    char *parc_data = (char *)hpx_parcel_get_data(p);
    hpx_addr_t *gmdata = reinterpret_cast<hpx_addr_t *>(parc_data);
    *gmdata = data;
    T **loc = reinterpret_cast<T **>(parc_data + sizeof(hpx_addr_t));
    *loc = location;

    // copy data into it
    void *arrdata = static_cast<void *>(parc_data + sizeof(T *)
                                          + sizeof(hpx_addr_t));
    for (size_t i = 0; i < local->local_count; ++i) {
      arrdata = local->manager->serialize(local->data + i, arrdata);
    }

    // send parcel
    hpx_parcel_set_action(p, array_collect_receive_);
    hpx_parcel_set_target(p, HPX_THERE(0));
    hpx_parcel_set_cont_action(p, hpx_lco_set_action);
    hpx_parcel_set_cont_target(p, done);
    hpx_parcel_send_sync(p);

    hpx_gas_unpin(global);

    return HPX_SUCCESS;
  }

  static int array_collect_receive_handler(char *data, size_t size) {
    hpx_addr_t *meta = reinterpret_cast<hpx_addr_t *>(data);
    hpx_addr_t global = hpx_addr_add(*meta,
              sizeof(ArrayMetaData<T>) * hpx_get_my_rank(),
              sizeof(ArrayMetaData<T>));
    ArrayMetaData<T> *local{nullptr};
    assert(hpx_gas_try_pin(global, (void **)&local));

    T *location = *(reinterpret_cast<T **>(data + sizeof(hpx_addr_t)));
    char *incoming = data + sizeof(hpx_addr_t) + sizeof(T *);
    char *final = data + size;
    size_t i{0};
    while (incoming != final) {
      incoming = (char *)local->manager->deserialize(incoming, location + i++);
    }

    hpx_gas_unpin(global);

    return HPX_SUCCESS;
  }

  static int array_set_manager_handler(Serializer *manager, hpx_addr_t data) {
    hpx_addr_t global = hpx_addr_add(data,
              sizeof(ArrayMetaData<T>) * hpx_get_my_rank(),
              sizeof(ArrayMetaData<T>));
    ArrayMetaData<T> *local{nullptr};
    assert(hpx_gas_try_pin(global, (void **)&local));
    if (local->manager != nullptr) {
      delete local->manager;
    }
    local->manager = manager;
    hpx_gas_unpin(global);
    hpx_exit(0, nullptr);
  }

};


template <typename T>
hpx_action_t Array<T>::array_local_count_ = HPX_ACTION_NULL;

template <typename T>
hpx_action_t Array<T>::array_total_count_ = HPX_ACTION_NULL;

template <typename T>
hpx_action_t Array<T>::allocate_array_meta_ = HPX_ACTION_NULL;

template <typename T>
hpx_action_t Array<T>::allocate_local_work_ = HPX_ACTION_NULL;

template <typename T>
hpx_action_t Array<T>::allocate_array_destroy_reducer_ = HPX_ACTION_NULL;

template <typename T>
hpx_action_t Array<T>::deallocate_array_ = HPX_ACTION_NULL;

template <typename T>
hpx_action_t Array<T>::deallocate_array_local_ = HPX_ACTION_NULL;

template <typename T>
hpx_action_t Array<T>::get_or_put_retcode_reducer_ = HPX_ACTION_NULL;

template <typename T>
hpx_action_t Array<T>::get_or_put_reducer_delete_ = HPX_ACTION_NULL;

template <typename T>
hpx_action_t Array<T>::array_get_ = HPX_ACTION_NULL;

template <typename T>
hpx_action_t Array<T>::array_put_ = HPX_ACTION_NULL;

template <typename T>
hpx_action_t Array<T>::array_segment_request_ = HPX_ACTION_NULL;

template <typename T>
hpx_action_t Array<T>::array_collect_ = HPX_ACTION_NULL;

template <typename T>
hpx_action_t Array<T>::array_collect_request_ = HPX_ACTION_NULL;

template <typename T>
hpx_action_t Array<T>::array_collect_receive_ = HPX_ACTION_NULL;

template <typename T>
hpx_action_t Array<T>::array_set_manager_ = HPX_ACTION_NULL;


} // namespace dashmm


#endif // __DASHMM_ARRAY_H__
