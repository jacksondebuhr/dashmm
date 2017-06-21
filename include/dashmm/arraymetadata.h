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


#ifndef __DASHMM_ARRAY_META_DATA_H__
#define __DASHMM_ARRAY_META_DATA_H__


/// \file
/// \brief Definition of ArrayMetaData.


#include <hpx/hpx.h>


namespace dashmm {


/// Meta data for the array object
///
/// Array objects in DASHMM are two part objects in the global address space.
/// The Array points to the meta data which is the information below.
/// Somewhere else in the global address space will be the data itself.
template <typename T>
struct ArrayMetaData {
  T *data;              /// address of local segment of the array
  size_t local_count;   /// the number of records in the local portion
  size_t total_count;   /// the total number of records in all portions
  size_t size;          /// the size (bytes) of each record
};


} // namespace dashmm


#endif // __DASHMM_ARRAY_META_DATA_H__
