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


uint32_t dashmm_object_class(dashmm_object_t *object);

bool dashmm_user_object(dashmm_object_t *object);

bool dashmm_system_object(dashmm_object_t *object);


#endif
