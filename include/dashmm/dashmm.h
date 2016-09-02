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


#ifndef __DASHMM_USER_INTERFACE_H__
#define __DASHMM_USER_INTERFACE_H__


// The basic interface
#include "dashmm/array.h"
#include "dashmm/evaluator.h"
#include "dashmm/initfini.h"
#include "dashmm/types.h"

// The built in methods
#include "builtins/bh_method.h"
#include "builtins/direct_method.h"
#include "builtins/fmm_method.h"
#include "builtins/fmm97_method.h"

// The built in expansions
#include "builtins/laplace_com.h"
#include "builtins/laplace_com_acc.h"
#include "builtins/laplace_sph.h"
#include "builtins/yukawa.h"

// The built in distribution policies
#include "builtins/singlelocdistro.h"
#include "builtins/randomdistro.h"


#endif // __DASHMM_USER_INTERFACE_H__
