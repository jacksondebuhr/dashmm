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

#ifndef __DASHMM_OBJECT_H__
#define __DASHMM_OBJECT_H__


typedef struct {
  uint32_t properties;
} dashmm_object_t;


#define DASHMM_CLASS_MASK    0xffff
#define DASHMM_PROPERTY_MASK 0xffff0000


//We define the object classes and properties here so that we can guarantee
// uniqueness of the identifiers.
#define DASHMM_CLASS_KERNEL  1
#define DASHMM_CLASS_METHOD  2
#define DASHMM_CLASS_ARRAY   3


#define DASHMM_PROPERTY_SYSTEM 0x10000
#define DASHMM_PROPERTY_USER   0x20000


bool dashmm_verify_user_object(dashmm_handle_t handle, uint32_t type);


#endif // __DASHMM_OBJECT_H__
