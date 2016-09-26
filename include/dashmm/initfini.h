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


#ifndef __DASHMM_INIT_FINI_H__
#define __DASHMM_INIT_FINI_H__


/// \file
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
/// \return kSuccess
ReturnCode finalize();


} // namespace dashmm


#endif // __DASHMM_INIT_FINI_H__
