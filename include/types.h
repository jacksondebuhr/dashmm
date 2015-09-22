#ifndef __DASHMM_TYPES_H__
#define __DASHMM_TYPES_H__


#include "hpx/hpx.h"


namespace dashmm {


//Used to give meaning to return values
enum ReturnCode {
  kSuccess = 0,
  kRuntimeError = 1,
  kIncompatible = 2,
  kAllocationError = 3,
  kInitError = 4,
  kFiniError = 5
  //etc...
};


using ObjectHandle = hpx_addr_t;


} // namespace dashmm


#endif // __DASHMM_TYPES_H__
