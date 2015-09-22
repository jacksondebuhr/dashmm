#ifndef __DASHMM_BASIC_INTERFACE_H__
#define __DASHMM_BASIC_INTERFACE_H__


#include "include/builtins.h"
#include "include/expansion.h"
#include "include/method.h"
#include "include/types.h"


namespace dashmm {


/// \brief initialize DASHMM
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
/// \return kSuccess on successful initialization; kInitError otherwise
int init(int *argc, char ***argv);


/// \breif finalize DASHMM
///
/// This will finalize the runtime system supporting DASHMM and will free any
/// resources claimed by DASHMM.
///
/// \return kSuccess on successful shutdown; kFiniError otherwise
int finalize();


//TODO;
int evaluate(ObjectHandle sources, ObjectHandle targets, int refinement_limit,
             Method *method, Expansion *expansion);


/// \brief Allocate an array in DASHMM's global address space
///
/// This will allocate an array of @param count records each of @param size
/// bytes. The object so created will be able to be referenced by the handle
/// @param obj.
///
/// \param count - the number of elements of the array
/// \param size - the size in bytes of each element
/// \param obj [out] - a handle that references the newly allocated array
///
/// \return kSuccess on success; kRuntimeError if there is an error with the
///         runtime; kAllocationError if there was a problem with the allocation
int allocate_array(size_t count, size_t size, ObjectHandle *obj);


/// \brief Frees an array in DASHMM's global address space
///
/// This will release the array referenced by @param obj. After this call,
/// further use of @param obj will cause undefined behavior.
///
/// \param obj - handle that references the array object
///
/// \return kSuccess on successful deallocation; kRuntimeError if there is an
///         error from the runtime
int deallocate_array(ObjectHandle obj);


/// \brief Put some data into a DASHMM array
///
/// This will copy the provided @param in_data into the array with handle
/// @param obj. The data will be copied into the range
/// [@param first, @param last) in the array.
///
/// \param obj - handle to the array
/// \param first - first index in the array to which the put will occur
/// \param last - one past the last index in the array to which the put will
///               occur
/// \param in_data - the data to put into the array
///
/// \return kSuccess on success; kRuntimeError otherwise
int array_put(ObjectHandle obj, size_t first, size_t last, void *in_data);


/// \brief Get some data from a DASHMM array
///
/// This will copy data from the array with handle @param obj into the buffer
/// specified by @param out_data. The data will be copied from the range
/// [@param first, @param last) in the array.
///
/// \param obj - handle to the array
/// \param first - first index in the array from which the get will occur
/// \param last - one past the last index in the array from which the put will
///               occur
/// \param in_data [out] - buffer for retrieved data
///
/// \return kSuccess on success; kRuntimeError otherwise
int array_get(ObjectHandle obj, size_t first, size_t last, void *out_data);


} // namespace dashmm


#endif // __DASHMM_BASIC_INTERFACE_H__
