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


#ifndef __DASHMM_USER_INTERFACE_H__
#define __DASHMM_USER_INTERFACE_H__


// The basic interface
#include "dashmm/array.h"
#include "dashmm/broadcast.h"
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
#include "builtins/laplace.h"
#include "builtins/yukawa.h"

// The built in distribution policies
#include "builtins/singlelocdistro.h"
#include "builtins/randomdistro.h"
#include "builtins/fmm97distro.h"

#endif // __DASHMM_USER_INTERFACE_H__
