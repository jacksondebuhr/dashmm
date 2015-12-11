#ifndef __DASHMM_ARRAY_H__
#define __DASHMM_ARRAY_H__


#include <cstdlib>


#include <hpx/hpx.h>


namespace dashmm {


struct ArrayMetaData {
  hpx_addr_t data;
  size_t count;
  size_t size;
};


} // namespace dashmm


#endif // __DASHMM_ARRAY_H__
