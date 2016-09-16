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


/// \file src/array.cc
/// \brief Implementation of DASHMM Array actions.


#include <cassert>
#include <cstring>

#include <hpx/hpx.h>

#include "dashmm/array.h"
#include "dashmm/reductionops.h"


namespace dashmm {

struct AddrMergeData {
  int count;
  hpx_addr_t addr[];
};

struct AddrMergeArgs {
  int rank;
  hpx_addr_t addx;
};


// Operations for the address sharing
void init_handler(AddrMergeData *head, size_t bytes,
                  AddrMergeData *init, size_t init_bytes) {
  head->count = (bytes - sizeof(AddrMergeData)) / sizeof(hpx_addr_t);
  for (int i = 0; i < head->count; ++i) {
    head->addr[i] = HPX_NULL;
  }
}
HPX_ACTION(HPX_FUNCTION, HPX_ATTR_NONE, init_action, init_handler,
           HPX_POINTER, HPX_SIZE_T, HPX_POINTER, HPX_SIZE_T);

void op_handler(AddrMergeData *lhs, AddrMergeArgs *rhs, size_t bytes) {
  assert(bytes == sizeof(AddrMergeArgs));
  assert(rhs->rank < hpx_get_num_ranks());
  lhs->addr[rhs->rank] = rhs->addx;
  lhs->count -= 1;
}
HPX_ACTION(HPX_FUNCTION, HPX_ATTR_NONE, op_action, op_handler,
           HPX_POINTER, HPX_POINTER, HPX_SIZE_T);

bool pred_handler(AddrMergeData *lhs, size_t bytes) {
  return lhs->count == 0;
}
HPX_ACTION(HPX_FUNCTION, HPX_ATTR_NONE, pred_action, pred_handler,
           HPX_POINTER, HPX_SIZE_T);


// allocate the array meta data
int allocate_array_meta_handler(void *UNUSED, size_t UNWANTED) {
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

  hpx_exit(sizeof(retval), &retval);
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
           allocate_array_meta_action, allocate_array_meta_handler,
           HPX_POINTER, HPX_SIZE_T);


// the local part of the allocation work
int allocate_local_work_handler(hpx_addr_t data, hpx_addr_t reducer,
                                size_t record_size, size_t record_count) {
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
  local->data = HPX_NULL;

  delete [] contrib;

  int retval{0};
  if (record_count) {
    local->data = hpx_gas_alloc_local(1, record_count * record_size, 0);
    if (local->data == HPX_NULL) {
      retval = kAllocationError;
    }
  }

  hpx_gas_unpin(global);
  hpx_exit(sizeof(retval), &retval);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           allocate_local_work_action, allocate_local_work_handler,
           HPX_ADDR, HPX_ADDR, HPX_SIZE_T, HPX_SIZE_T);


// delete the reducer
int allocate_array_destroy_reducer_handler(hpx_addr_t reducer) {
  hpx_lco_delete_sync(reducer);
  hpx_exit(0, nullptr);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           allocate_array_destroy_reducer_action,
           allocate_array_destroy_reducer_handler,
           HPX_ADDR);


// Delete the local portion of the array
int deallocate_array_local_handler(hpx_addr_t meta) {
  hpx_addr_t global = hpx_addr_add(meta,
                                   sizeof(ArrayMetaData) * hpx_get_my_rank(),
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  assert(hpx_gas_try_pin(global, (void **)&local));

  if (local->data != HPX_NULL) {
    hpx_gas_free_sync(local->data);
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
HPX_ACTION(HPX_DEFAULT, 0, deallocate_array_action, deallocate_array_handler,
           HPX_ADDR);


/// Put data into an array object
///
/// This will fill records in a global array with the provided data. This
/// action is not well behaved as regards bad input arguments. This will likely
/// be updated in the future.
///
/// \param obj - the global address of the ArrayMetaData
/// \param first - the first record to put data into
/// \param last - the last (exclusive) record to put data into
/// \param in_data - the data to copy into the records
///
/// \returns - HPX_SUCCESS
int array_put_handler(hpx_addr_t obj, size_t first, size_t last,
                      void *in_data) {
  int retval = kSuccess;

  hpx_addr_t global = hpx_addr_add(obj,
                                   sizeof(ArrayMetaData) * hpx_get_my_rank(),
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  if (!hpx_gas_try_pin(global, (void **)&local)) {
    retval = kRuntimeError;
    hpx_exit(sizeof(retval), &retval);
  }

  char *data{nullptr};
  if (local->data == HPX_NULL) {
    hpx_gas_unpin(global);
    hpx_exit(sizeof(retval), &retval);
  }
  if (!hpx_gas_try_pin(local->data, (void **)&data)) {
    hpx_gas_unpin(global);
    retval = kRuntimeError;
    hpx_exit(sizeof(retval), &retval);
  }

  // Do some simple bounds checking
  if (last > local->local_count || last < first) {
    retval = kDomainError;
    hpx_gas_unpin(local->data);
    hpx_gas_unpin(global);
    hpx_exit(sizeof(retval), &retval);
  }

  char *beginning = data + first * local->size;
  size_t copy_size = (last - first) * local->size;
  memcpy(beginning, in_data, copy_size);

  hpx_gas_unpin(local->data);
  hpx_gas_unpin(global);

  hpx_exit(sizeof(retval), &retval);
}
HPX_ACTION(HPX_DEFAULT, 0, array_put_action, array_put_handler,
           HPX_ADDR, HPX_SIZE_T, HPX_SIZE_T, HPX_POINTER);


/// Get data from an array object
///
/// This will read records from a global array into the provided buffer. This
/// action is not well behaved as regards bad input arguments. This will likely
/// be updated in the future.
///
/// \param obj - the global address of the ArrayMetaData
/// \param first - the first record to read
/// \param last - the last (exclusive) record to read
/// \param out_data - a buffer into which the read data is stored
///
/// \returns - HPX_SUCCESS
int array_get_handler(hpx_addr_t obj, size_t first, size_t last,
                      void *out_data) {
  int retval = kSuccess;

  hpx_addr_t global = hpx_addr_add(obj,
                                   sizeof(ArrayMetaData) * hpx_get_my_rank(),
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  if (!hpx_gas_try_pin(global, (void **)&local)) {
    retval = kRuntimeError;
    hpx_exit(sizeof(retval), &retval);;
  }

  char *data{nullptr};
  if (local->data == HPX_NULL) {
    hpx_gas_unpin(global);
    hpx_exit(sizeof(retval), &retval);;
  }
  if (!hpx_gas_try_pin(local->data, (void **)&data)) {
    retval = kRuntimeError;
    hpx_gas_unpin(global);
    hpx_exit(sizeof(retval), &retval);;
  }

  // Do some simple bounds checking
  if (last > local->local_count || last < first) {
    retval = kDomainError;
    hpx_gas_unpin(local->data);
    hpx_gas_unpin(global);
    hpx_exit(sizeof(retval), &retval);;
  }

  char *beginning = data + first * local->size;
  size_t copy_size = (last - first) * local->size;
  memcpy(out_data, beginning, copy_size);

  hpx_gas_unpin(local->data);
  hpx_gas_unpin(global);

  hpx_exit(sizeof(retval), &retval);;
}
HPX_ACTION(HPX_DEFAULT, 0, array_get_action, array_get_handler,
           HPX_ADDR, HPX_SIZE_T, HPX_SIZE_T, HPX_POINTER);


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


int array_collect_prep_handler(hpx_addr_t UNUSED) {
  // This is used to signal that prep is done
  hpx_addr_t retval[2];
  int n_ranks = hpx_get_num_ranks();
  retval[0] = hpx_lco_reduce_new(n_ranks, sizeof(size_t) * n_ranks,
                                 size_sum_ident, size_sum_op);
  int nonsense{0};
  retval[1] = hpx_lco_user_new(
        sizeof(hpx_addr_t) * n_ranks + sizeof(AddrMergeData), init_action,
        op_action, pred_action, &nonsense, sizeof(nonsense));
  hpx_exit(sizeof(hpx_addr_t) * 2, &retval);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           array_collect_prep_action, array_collect_prep_handler,
           HPX_ADDR);


int array_collect_handler(hpx_addr_t data, hpx_addr_t offsets,
                          hpx_addr_t addxes) {
  // get local meta data
  int my_rank = hpx_get_my_rank();
  int n_ranks = hpx_get_num_ranks();

  hpx_addr_t global = hpx_addr_add(data,
                                   sizeof(ArrayMetaData) * my_rank,
                                   sizeof(ArrayMetaData));
  ArrayMetaData *local{nullptr};
  assert(hpx_gas_try_pin(global, (void **)&local));

  // Reduce the offsets
  size_t *contrib = new size_t[n_ranks];
  for (int i = 0; i < n_ranks; ++i) {
    if (i >= my_rank) {
      contrib[i] = local->local_count;
    } else {
      contrib[i] = 0;
    }
  }
  hpx_lco_set_lsync(offsets, sizeof(size_t) * n_ranks, contrib, HPX_NULL);

  // reduce the data locations
  AddrMergeArgs args{my_rank, local->data};
  hpx_lco_set_lsync(addxes, sizeof(args), &args, HPX_NULL);

  // Rank zero will then collect the information
  char *retval{nullptr};
  if (my_rank == 0) {
    retval = new char[local->size * local->total_count];

    // compute offsets
    hpx_lco_get(offsets, sizeof(size_t) * n_ranks, contrib);
    size_t *counts = new size_t[n_ranks];
    counts[0] = contrib[0];
    for (int i = 1; i < n_ranks; ++i) {
      counts[i] = contrib[i] - contrib[i - 1];
    }

    // loop over ranks memget
    size_t datasize = sizeof(AddrMergeData) + sizeof(hpx_addr_t) * n_ranks;
    AddrMergeData *where =
        reinterpret_cast<AddrMergeData *>(new char[datasize]);
    hpx_lco_get(addxes, datasize, where);

    hpx_gas_memget_sync(retval, where->addr[0], local->size * counts[0]);
    for (int i = 1; i < n_ranks; ++i) {
      hpx_gas_memget_sync(&retval[contrib[i-1] * local->size],
                          where->addr[i],
                          local->size * counts[i]);
    }

    delete [] counts;
    delete [] where;
  }

  delete [] contrib;
  hpx_gas_unpin(global);
  hpx_exit(sizeof(retval), &retval);
}
HPX_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
           array_collect_action, array_collect_handler,
           HPX_ADDR, HPX_ADDR, HPX_ADDR);


} // namespace dashmm
