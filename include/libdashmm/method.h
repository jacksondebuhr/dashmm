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

#ifndef __DASHMM_METHOD_H__
#define __DASHMM_METHOD_H__


#include "libdashmm/basic_types.h"
#include "libdashmm/object.h"


///
/// \type dashmm_method_t
///
/// This type specifies the details of the method to use. This defines not
/// only the overall class of method (FMM or BH), but details about the 
/// method, such as the cell acceptance criteria for the given method.
///
/// JD: This will not be exposed. The object will be dealt with through an
///     interface.
///
/// JD: We need to define what sort of functions are used for the MAC. We should
///     give these as a typedef.
///
typedef struct {
  dashmm_object_t object;
  
  int use_local_operators;
  int use_MAC;
  hpx_action_t MAC;
  //These are planned future extensions and will not be useful for anything
  // immediately.
  int use_merge_and_shift;
  int use_irregular_trees;
} dashmm_method_t;


#endif
