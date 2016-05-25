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


#ifndef __DASHMM_INIT_FINI_H__
#define __DASHMM_INIT_FINI_H__


/// \file include/dashmm/basic.h
/// \brief The basic interface to DASHMM.


#include "dashmm/types.h"


/// \namespace dashmm
/// \brief The namespace inside which all DASHMM symbols are defined
namespace dashmm {


/// initialize DASHMM
///
/// This will initialize the runtime system supporting DASHMM and will setup
/// any resources that DASHMM requires. The runtime system can have some of
/// its behavior modified by command line arguments. Any arguments that the
/// runtime uses will be removed from the list of arguments. After this call
/// the remaining arguments (which should be application specific) can be
/// processed by the application.
///
/// \param argc [inout] - the number of command line arguments
/// \param argv [inout] - the arguments themselves
///
/// \return kSuccess on successful initialization; kRuntimeError if there is a
///         problem initializing the HPX-5 runtime;
ReturnCode init(int *argc, char ***argv);


/// finalize DASHMM
///
/// This will finalize the runtime system supporting DASHMM and will free any
/// resources claimed by DASHMM.
///
/// \return kSuccess on successful shutdown; kFiniError otherwise
ReturnCode finalize();


} // namespace dashmm


#endif // __DASHMM_INIT_FINI_H__
