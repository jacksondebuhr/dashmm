#ifndef __DASHMM_TYPES_H__
#define __DASHMM_TYPES_H__


/// \file include/types.h
/// \brief Basic type definitions used throughout DASHMM


#include "hpx/hpx.h"


namespace dashmm {


/// Return codes from DASHMM library calls
///
/// The possible return codes from DASHMM calls. For details about the
/// meanings of these, please see the individual library calls that generate
/// them.
enum ReturnCode {
  kSuccess = 0,
  kRuntimeError = 1,
  kIncompatible = 2,
  kAllocationError = 3,
  kInitError = 4,
  kFiniError = 5
  //etc...
};


/// \class ObjectHandle
/// \brief An opaque handle to global objects in DASHMM
///
/// References to objects created by DASHMM that are returned to the user
/// will all be of class ObjectHandle. These references are opaque, and there
/// are no methods exposed to the user.
using ObjectHandle = hpx_addr_t;


} // namespace dashmm


#endif // __DASHMM_TYPES_H__
