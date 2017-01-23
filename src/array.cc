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


/// \file
/// \brief Implementation of DASHMM Array actions.


#include <cassert>
#include <cstring>

#include <hpx/hpx.h>

#include "dashmm/array.h"
#include "dashmm/reductionops.h"


namespace dashmm {


/// Data stored in the address merge reduction LCO
struct AddrMergeData {
  int count;
  hpx_addr_t addr[];
};

/// Data required for set in the address merge reduction LCO
struct AddrMergeArgs {
  int rank;
  hpx_addr_t addx;
};


/// Initialization operation for address merge reduction LCO
void init_handler(AddrMergeData *head, size_t bytes,
                  AddrMergeData *init, size_t init_bytes) {
  head->count = (bytes - sizeof(AddrMergeData)) / sizeof(hpx_addr_t);
  for (int i = 0; i < head->count; ++i) {
    head->addr[i] = HPX_NULL;
  }
}
HPX_ACTION(HPX_FUNCTION, HPX_ATTR_NONE, init_action, init_handler,
           HPX_POINTER, HPX_SIZE_T, HPX_POINTER, HPX_SIZE_T);

/// Reduction operation for address merge reduction LCO
void op_handler(AddrMergeData *lhs, AddrMergeArgs *rhs, size_t bytes) {
  assert(bytes == sizeof(AddrMergeArgs));
  assert(rhs->rank < hpx_get_num_ranks());
  lhs->addr[rhs->rank] = rhs->addx;
  lhs->count -= 1;
}
HPX_ACTION(HPX_FUNCTION, HPX_ATTR_NONE, op_action, op_handler,
           HPX_POINTER, HPX_POINTER, HPX_SIZE_T);

/// Predicate operation for address merge reduction LCO
bool pred_handler(AddrMergeData *lhs, size_t bytes) {
  return lhs->count == 0;
}
HPX_ACTION(HPX_FUNCTION, HPX_ATTR_NONE, pred_action, pred_handler,
           HPX_POINTER, HPX_SIZE_T);


/// Action for allocating Array metadata
///
/// This action is called on a single rank to allocate some collective
/// data that is needed to allocate an array. The returned information is
/// broadcast back to every rank during hpx_exit().
///
/// \returns - HPX_SUCCESS
int allocate_array_meta_handler() {
  ArrayMetaAllocRunReturn retval{HPX_NULL, HPX_NULL, kSuccess};

  int ranks = hpx_get_num_ranks();
  retval.meta = hpx_gas_alloc_cyclic(ranks, sizeof(ArrayMetaData), 0);
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
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           allocate_array_meta_action, allocate_array_meta_handler);


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
int allocate_local_work_handler(hpx_addr_t data, hpx_addr_t reducer,
                                hpx_addr_t retcode,
                                size_t record_size, size_t record_count,
                                char *segment) {
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
                                   hpx_get_my_rank() * sizeof(ArrayMetaData),
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  assert(hpx_gas_try_pin(global, (void **)&local));
  local->local_count = record_count;
  local->total_count = contrib[ranks - 1];
  local->size = record_size;
  local->data = segment;

  delete [] contrib;

  int retval{0};
  if (record_count && local->data == nullptr) {
    try {
      local->data = new char [record_count * record_size];
    } catch (std::bad_alloc &ba) {
      retval = kAllocationError;
    }
  }

  hpx_lco_set_lsync(retcode, sizeof(int), &retval, HPX_NULL);
  hpx_lco_get(retcode, sizeof(int), &retval);

  hpx_gas_unpin(global);
  hpx_exit(sizeof(retval), &retval);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           allocate_local_work_action, allocate_local_work_handler,
           HPX_ADDR, HPX_ADDR, HPX_ADDR, HPX_SIZE_T, HPX_SIZE_T,
           HPX_POINTER);


/// Action to deallocate the reduction LCO used during Array allocation
///
/// \param reducer - the global address of the reduction LCO
///
/// \returns - HPX_SUCCESS
int allocate_array_destroy_reducer_handler(hpx_addr_t reducer,
                                           hpx_addr_t retcode) {
  hpx_lco_delete_sync(reducer);
  hpx_lco_delete_sync(retcode);
  hpx_exit(0, nullptr);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           allocate_array_destroy_reducer_action,
           allocate_array_destroy_reducer_handler,
           HPX_ADDR, HPX_ADDR);


/// Action that deletes the local portion of an Array
///
/// This action is the target of a broadcast. It will delete the local portion
/// of the Array's global memory.
///
/// \param meta - the global address of the Array's meta data.
///
/// \returns - HPX_SUCCESS
int deallocate_array_local_handler(hpx_addr_t meta) {
  hpx_addr_t global = hpx_addr_add(meta,
                                   sizeof(ArrayMetaData) * hpx_get_my_rank(),
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  assert(hpx_gas_try_pin(global, (void **)&local));

  if (local->data != HPX_NULL) {
    delete [] local->data;
  }

  hpx_gas_unpin(global);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           deallocate_array_local_action, deallocate_array_local_handler,
           HPX_ADDR);


/// Action that deallocates an array object
///
/// This will delete the ArrayMetaData and the records themselves.
///
/// \param obj - the global address of the array's meta data
///
/// \returns - HPX_SUCCESS
int deallocate_array_handler(hpx_addr_t obj) {
  hpx_bcast_rsync(deallocate_array_local_action, &obj);
  hpx_gas_free_sync(obj);
  hpx_exit(0, nullptr);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           deallocate_array_action, deallocate_array_handler,
           HPX_ADDR);


/// Action that allocates a reducer for return codes
///
/// To be very careful, we collect the returns codes from each rank and
/// take the maximum value.
///
/// \returns - HPX_SUCCESS
int get_or_put_retcode_reducer_handler(void) {
  hpx_addr_t retval = hpx_lco_reduce_new(hpx_get_num_ranks(),
      sizeof(int), int_max_ident_op, int_max_op);
  hpx_exit(sizeof(retval), &retval);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           get_or_put_retcode_reducer_action,
           get_or_put_retcode_reducer_handler);


/// Action that deletes a reducer for return codes
///
/// \param lco - the LCO's address
///
/// \returns - HPX_SUCCESS
int get_or_put_reducer_delete_handler(hpx_addr_t lco) {
  hpx_lco_delete_sync(lco);
  hpx_exit(0, nullptr);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           get_or_put_reducer_delete_action,
           get_or_put_reducer_delete_handler,
           HPX_ADDR);


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
int array_put_handler(hpx_addr_t obj, size_t first, size_t last,
                      void *in_data, hpx_addr_t reducer) {
  int retval = kSuccess;

  hpx_addr_t global = hpx_addr_add(obj,
                                   sizeof(ArrayMetaData) * hpx_get_my_rank(),
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  if (!hpx_gas_try_pin(global, (void **)&local)) {
    retval = kRuntimeError;
  } else {
    if (last > local->local_count || last < first) {
      retval = kDomainError;
    } else if (local->local_count) {
      assert(local->data != nullptr);

      char *beginning = local->data + first * local->size;
      size_t copy_size = (last - first) * local->size;
      memcpy(beginning, in_data, copy_size);
    }

    hpx_gas_unpin(global);
  }

  hpx_lco_set_lsync(reducer, sizeof(int), &retval, HPX_NULL);
  hpx_lco_get(reducer, sizeof(int), &retval);

  hpx_exit(sizeof(retval), &retval);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE, array_put_action, array_put_handler,
           HPX_ADDR, HPX_SIZE_T, HPX_SIZE_T, HPX_POINTER, HPX_ADDR);


/// Get data from an array object
///
/// This will read records from a global array into @p out_data. This
/// action is not well behaved as regards bad input arguments. This will likely
/// be updated in the future.
///
/// \param obj - the global address of the ArrayMetaData
/// \param first - the first record to read
/// \param last - the last (exclusive) record to read
/// \param out_data - a buffer into which the read data is stored
/// \param reducer - the LCO that will reduce the error codes from each rank
///
/// \returns - HPX_SUCCESS
int array_get_handler(hpx_addr_t obj, size_t first, size_t last,
                      void *out_data, hpx_addr_t reducer) {
  int retval = kSuccess;

  hpx_addr_t global = hpx_addr_add(obj,
                                   sizeof(ArrayMetaData) * hpx_get_my_rank(),
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  if (!hpx_gas_try_pin(global, (void **)&local)) {
    retval = kRuntimeError;
  } else {
    if (last > local->local_count || last < first) {
      retval = kDomainError;
    } else if (local->local_count) {
      assert(local->data != nullptr);

      char *beginning = local->data + first * local->size;
      size_t copy_size = (last - first) * local->size;
      memcpy(out_data, beginning, copy_size);
    }

    hpx_gas_unpin(global);
  }

  // reduce on error condition
  hpx_lco_set_lsync(reducer, sizeof(int), &retval, HPX_NULL);
  hpx_lco_get(reducer, sizeof(int), &retval);

  hpx_exit(sizeof(retval), &retval);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE, array_get_action, array_get_handler,
           HPX_ADDR, HPX_SIZE_T, HPX_SIZE_T, HPX_POINTER, HPX_ADDR);


/// Action to return the local length of the Array
///
/// \param data - global address of the Array's meta data
///
/// \returns - HPX_SUCCESS
int array_local_count_handler(hpx_addr_t data) {
  hpx_addr_t global = hpx_addr_add(data,
                                   sizeof(ArrayMetaData) * hpx_get_my_rank(),
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  assert(hpx_gas_try_pin(global, (void **)&local));

  size_t retval = local->local_count;

  hpx_gas_unpin(global);
  hpx_exit(sizeof(retval), &retval);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           array_local_count_action, array_local_count_handler,
           HPX_ADDR);


/// Action to return the Array's total length
///
/// \param data - global address of the Array's meta data
///
/// \returns - HPX_SUCCESS
int array_total_count_handler(hpx_addr_t data) {
  hpx_addr_t global = hpx_addr_add(data,
                                   sizeof(ArrayMetaData) * hpx_get_my_rank(),
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  assert(hpx_gas_try_pin(global, (void **)&local));

  size_t retval = local->total_count;

  hpx_gas_unpin(global);
  hpx_exit(sizeof(retval), &retval);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           array_total_count_action, array_total_count_handler,
           HPX_ADDR);


static int array_collect_receive_handler(char *data, size_t size) {
  char *location = *(reinterpret_cast<char **>(data));
  memcpy(location, data + sizeof(char *), size - sizeof(char *));
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
           array_collect_receive_action,
           array_collect_receive_handler,
           HPX_POINTER, HPX_SIZE_T);


static int array_collect_request_handler(hpx_addr_t data, hpx_addr_t done,
                                         char *location) {
  hpx_addr_t global = hpx_addr_add(data,
                                   sizeof(ArrayMetaData) * hpx_get_my_rank(),
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  assert(hpx_gas_try_pin(global, (void **)&local));

  // get a parcel of the right size
  size_t arrsize = local->size * local->local_count;
  size_t msgsize = sizeof(char *) + arrsize;
  hpx_parcel_t *p = hpx_parcel_acquire(nullptr, msgsize);
  char *parc_data = (char *)hpx_parcel_get_data(p);
  char **loc = reinterpret_cast<char **>(parc_data);
  *loc = location;

  // copy data into it
  char *arrdata = parc_data + sizeof(char *);
  memcpy(arrdata, local->data, arrsize);

  // send parcel
  hpx_parcel_set_action(p, array_collect_receive_action);
  hpx_parcel_set_target(p, HPX_THERE(0));
  hpx_parcel_set_cont_action(p, hpx_lco_set_action);
  hpx_parcel_set_cont_target(p, done);
  hpx_parcel_send_sync(p);

  hpx_gas_unpin(global);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           array_collect_request_action,
           array_collect_request_handler,
           HPX_ADDR, HPX_ADDR, HPX_POINTER);


int array_collect_handler(hpx_addr_t data) {
  // collect counts from GAS
  int my_rank = hpx_get_my_rank();
  int n_rank = hpx_get_num_ranks();

  ArrayMetaData *meta = (ArrayMetaData *)hpx_malloc_registered(
                                            sizeof(ArrayMetaData) * n_rank);
  for (int i = 0; i < n_rank; ++i) {
    hpx_addr_t target = hpx_addr_add(data, sizeof(ArrayMetaData) * i,
                                     sizeof(ArrayMetaData));
    hpx_gas_memget_sync(&meta[i], target, sizeof(ArrayMetaData));
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
  char *retval{nullptr};
  retval = new char[meta[my_rank].size * meta[my_rank].total_count];

  // call on each rank with location and and gate and data
  for (int i = 0; i < n_rank; ++i) {
    if (i != my_rank) {
      char *location = retval + meta[i].size * offsets[i];
      hpx_call(HPX_THERE(i), array_collect_request_action, HPX_NULL,
               &data, &done, &location);
    }
  }

  memcpy(retval + meta[my_rank].size * offsets[my_rank],
         meta[my_rank].data,
         meta[my_rank].size * meta[my_rank].local_count);
  hpx_lco_and_set(done, HPX_NULL);

  // wait for results
  hpx_lco_wait(done);
  hpx_lco_delete_sync(done);

  hpx_free_registered(meta);
  delete [] offsets;

  hpx_exit(sizeof(retval), &retval);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           array_collect_action, array_collect_handler,
           HPX_ADDR);


int array_segment_request_handler(hpx_addr_t data) {
  hpx_addr_t global = hpx_addr_add(data,
                                   sizeof(ArrayMetaData) * hpx_get_my_rank(),
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  assert(hpx_gas_try_pin(global, (void **)&local));

  // collect info
  SegmentReturn retval{};
  retval.segment = local->data;
  retval.count = local->local_count;

  hpx_gas_unpin(global);

  hpx_exit(sizeof(retval), &retval);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           array_segment_request_action, array_segment_request_handler,
           HPX_ADDR);


} // namespace dashmm
