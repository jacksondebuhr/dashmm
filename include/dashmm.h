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

#ifndef __DASHMM_INTERFACE_H__
#define __DASHMM_INTERFACE_H__


#include "libdashmm/basic_types.h"
//TODO: a list of what types and so on are defined by this inclusion, with a 
//      note that we are purposefully hiding some of the runtime details from
//      those reading this file.
// dashmm_status_t
// dashmm_result_options_t
// dashmm_handle_t; INVALID_HANDLE


// ==================================================================
// Built-in methods and kernels
// ==================================================================

///
/// \var DASHMM_BARNES_HUT_METHOD
///
/// This is the handle for the built in Barnes-Hut method. This is as close
/// to the classical Barnes-Hut method as is available in DASHMM.
///
extern dashmm_handle_t DASHMM_BARNES_HUT_METHOD;

///
/// \var DASHMM_FAST_MULTIPOLE_METHOD
///
/// This is the handle for the built in FMM. This is as close to the classical
/// FMM that is available in DASHMM.
///
extern dashmm_handle_t DASHMM_FAST_MULTIPOLE_METHOD;

///
/// \var DASHMM_LAPLACE_POTENTIAL
///
/// This is the handle to the built-in Laplace kernel.
///
extern dashmm_handle_t DASHMM_LAPLACE_POTENTIAL;

///
/// \var DASHMM_YUKAWA_POTENTIAL
///
/// This is the handle to the build in Yukawa kernel.
///
extern dashmm_handle_t DASHMM_YUKAWA_POTENTIAL;


// ==================================================================
// Basic Interface to DASHMM
// ==================================================================



// ==================================================================
// Advanced Interface to DASHMM
// ==================================================================
//      The intent is that the user will need only to include a single file
//      to use the library. As such, the full interface will appear in this
//      header.


#endif
