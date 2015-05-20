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

#ifndef __DASHMM_KERNEL_H__
#define __DASHMM_KERNEL_H__


#include "libdashmm/basic_types.h"
#include "libdashmm/object.h"


///
/// \type dashmm_kernel_t
///
/// This type specifies the kernel information. This includes not only the 
/// particular potential, and the various translation operators, but also the
/// parameters in the potential.
///
/// JD: This will not be exposed. The object will be dealt with through an
///     interface.
///
/// JD: At the moment, the plan is that each kernel is complete for a given
///     problem and a given accuracy. So the Laplace kernel is actually a set
///     of kernels for each accuracy. These may be largely similar to one 
///     another, but they may not be.
///
/// JD: We need to define what sort of arguments these operators need, and
///     define these things.
///
typedef struct {
  dashmm_object_t object;
  
  hpx_action_t source_to_multiple;
  hpx_action_t multipole_to_multipole;
  hpx_action_t multipole_to_target;
  hpx_action_t multipole_to_local;
  hpx_action_t local_to_local;
  hpx_action_t local_to_target;
  //more for the merge-and-shift
  
  //various tables for maps into different accuracy levels?
  
  //perhaps also a function for creating the coefficients that are needed
  // by the method.
  
  //A list of values that just have to be provided at runtime? Or rather the
  // number of these. These would be constant over the operator uses during 
  // the evaluation. 
} dashmm_kernel_t;


#endif
