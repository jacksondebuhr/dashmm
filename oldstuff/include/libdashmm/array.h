// =============================================================================
//  DASHMM
//
//  Copyright (c) 2014 - 2015, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license.  See the LICENSE file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
//
//  Authors:
//    Jackson DeBuhr, Indiana University <jdebuhr [at] indiana.edu>
// =============================================================================

#ifndef __DASHMM_ARRAY_H__
#define __DASHMM_ARRAY_H__


#include "libdashmm/basic_types.h"
#include "libdashmm/object.h"


typedef enum {
  DASHMM_DISTRIB_CYCLIC,
  DASHMM_DISTRIB_BLOCK_CYCLIC,
  DASHMM_DISTRIB_BLOCKED,
} dashmm_array_distrib_t;


typedef struct {
  dashmm_object_t object;
  
  hpx_addr_t data;          // hpx address of the data
  size_t num_records;       // the number of records total
  size_t record_size;       // the size of each record
  size_t block_size;        // the number of records per block in GAS
} dashmm_array_t;


#endif //__DASHMM_ARRAY_H__
