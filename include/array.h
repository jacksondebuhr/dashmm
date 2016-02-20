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


extern hpx_action_t allocate_array_action;
extern hpx_action_t deallocate_array_action;
extern hpx_action_t array_put_action;
extern hpx_action_t array_get_action;



template <typename T>
class Array {
 public:
  Array(size_t record_count) {
    // TODO: Make work also from HPX thread.
    int runcode;
    int *arg = &runcode;
    size_t size = sizeof(T);
    hpx_run(&allocate_array_action, &record_count, &size, &data_, &arg);
  }

  Array(hpx_addt_t data) : data_{data} { }

  bool valid() const {return data != HPX_NULL;}

  hpx_addr_t data() const {return data_;}

  // TODO: Do we want queries on the number of records?

  ReturnCode destroy() {
    // TODO: Make work also from HPX thread.
    if (HPX_SUCCESS == hpx_run(&deallocate_array_action, &data_)) {
      return kSuccess;
    } else {
      return kRuntimeError;
    }
  };

  ReturnCode get(size_t first, size_t last, T *out_data) {
    // TODO: Make work also from HPX thread.
    int runcode;
    int *arg = &runcode;
    hpx_run(&array_get_action, &data_, &first, &last, &arg, &out_data);
    return static_cast<ReturnCode>(runcode);
  }

  ReturnCode put(size_t first, size_t last, T *in_data) {
    // TODO: Make work also from HPX thread.
    int runcode;
    int *arg = &runcode;
    hpx_run(&array_put_action, &data_, &first, &last, &arg, &in_data);
    return static_cast<ReturnCode>(runcode);
  }

  // TODO: do we define an iterator for this? What does that even mean in the
  // context of GAS?
 private:
  hpx_addr_t data_;
};



} // namespace dashmm


#endif // __DASHMM_ARRAY_H__
