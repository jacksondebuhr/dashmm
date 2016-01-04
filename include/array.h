#ifndef __DASHMM_ARRAY_H__
#define __DASHMM_ARRAY_H__


/// \file include/array.h
/// \brief Definitions needed to interact with DASHMM array objects.


#include <cstdlib>


#include <hpx/hpx.h>


namespace dashmm {


/// Meta data for the array object
///
/// Array objects in DASHMM are two part objects in the global address space.
/// The ObjectHandle points to the meta data which is the information below.
/// Somewhere else in the global address space will be the data itself.
/// Currently, as DASHMM is specialized to SMP operation, the array data will
/// be on the same locality. In the future, this may change, and so we separate
/// the array meta data from the array itself, as the required metadata may
/// need to change as the data layout becomes more flexible.
struct ArrayMetaData {
  hpx_addr_t data;    /// global address of the array data
  size_t count;       /// the number of records in the array
  size_t size;        /// the size (bytes) of each record
};


} // namespace dashmm


#endif // __DASHMM_ARRAY_H__
