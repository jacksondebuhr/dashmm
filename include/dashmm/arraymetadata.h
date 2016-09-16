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


#ifndef __DASHMM_ARRAY_META_DATA_H__
#define __DASHMM_ARRAY_META_DATA_H__


/// \file include/dashmm/arraymetadata.h
/// \brief Definition of ArrayMetaData.


#include <hpx/hpx.h>


namespace dashmm {


/// Meta data for the array object
///
/// Array objects in DASHMM are two part objects in the global address space.
/// The ObjectHandle points to the meta data which is the information below.
/// Somewhere else in the global address space will be the data itself.
struct ArrayMetaData {
  hpx_addr_t data;      /// global address of the array data
  size_t local_count;   /// the number of records in the local portion
  size_t total_count;   /// the total number of records in all portions
  size_t size;          /// the size (bytes) of each record
};


} // namespace dashmm


#endif // __DASHMM_ARRAY_META_DATA_H__
