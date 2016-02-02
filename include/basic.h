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


#ifndef __DASHMM_BASIC_INTERFACE_H__
#define __DASHMM_BASIC_INTERFACE_H__


/// \file include/basic.h
/// \brief The basic interface to DASHMM.


#include <memory>

#include "include/expansion.h"
#include "include/method.h"
#include "include/types.h"


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
///         problem initializing the HPX-5 runtime; kInitError otherwise
ReturnCode init(int *argc, char ***argv);


/// finalize DASHMM
///
/// This will finalize the runtime system supporting DASHMM and will free any
/// resources claimed by DASHMM.
///
/// \return kSuccess on successful shutdown; kFiniError otherwise
ReturnCode finalize();


/// Perform a multipole method evaluation
///
/// This is the central routine in basic DASHMM use. The user supplies a
/// @param method and an @param expansion along with source and target
/// information and DASHMM will compute the potential for the given
/// situation.
///
/// The source description is provided by giving an object handle to a
/// DASHMM array object as well as the offsets to the position and the charge
/// in each record of the sources. The positions are assumed to be three
/// double values and the charge is assumed to be a single double value.
///
/// The target description is provided by giving the object handle to the
/// DASHMM array holding the target data as well as the offset in each record
/// of that array to the position and the potential. The position is taken to
/// be three double values and the potential will be returned as a pair of
/// double values giving both the real and imaginary component of the complex
/// potential. For expansions that return only real potentials, the second
/// double in phi is not used, but the system will use it (to fill the value
/// with zero).
///
/// DASHMM will not modify the source array except in the case where the source
/// and target array are the same, and in that case only the potential field
/// will be modified.
///
/// The @param refinement_limit tells DASHMM at what point to stop partitioning
/// the source and target trees. If there are more than the given limit sources
/// or targets in a node, that node will be subdivided.
///
/// \param sources - handle to the array object holding the source data
/// \param spos_offset - offset in the records in the source array to the
///                      position of the sources.
/// \param q_offset - offset in the records in the source array to the charge of
///                   of the source
/// \param targets - handle to the array object holding the target data
/// \param tpos_offset - offset in the target record to the position of the
///                      target
/// \param phi_offset - offset in the target record to the potential for the
///                     target
/// \param refinement_limit - the limiting number of sources or targets in each
///                     leaf of the hierarchical subdivision of the source or
///                     target points
/// \param method - the method to use for the evaluation
/// \param expansion - the expansion to use for the evaluation
///
/// \returns - kSuccess on successful evaluation; kIncompatible if the given
///            method and expansion are incompatible; kRuntimeError if there
///            is an error from the runtime.
ReturnCode evaluate(ObjectHandle sources, size_t spos_offset, size_t q_offset,
                    ObjectHandle targets, size_t tpos_offset, size_t phi_offset,
                    int refinement_limit,
                    std::unique_ptr<Method> method,
                    std::unique_ptr<Expansion> expansion);


/// Allocate an array in DASHMM's global address space
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
ReturnCode allocate_array(size_t count, size_t size, ObjectHandle *obj);


/// Frees an array in DASHMM's global address space
///
/// This will release the array referenced by @param obj. After this call,
/// further use of @param obj will cause undefined behavior.
///
/// \param obj - handle that references the array object
///
/// \return kSuccess on successful deallocation; kRuntimeError if there is an
///         error from the runtime
ReturnCode deallocate_array(ObjectHandle obj);


/// Put some data into a DASHMM array
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
/// \return kSuccess on success; kRuntimeError if there are errors from the
///         runtime or kDomainError if the input arguments are invalid.
ReturnCode array_put(ObjectHandle obj, size_t first, size_t last,
                     void *in_data);


/// Get some data from a DASHMM array
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
/// \return kSuccess on success; kDomainError if the input arguments are invalid
///         of kRuntimeError otherwise.
ReturnCode array_get(ObjectHandle obj, size_t first, size_t last,
                     void *out_data);


} // namespace dashmm


#endif // __DASHMM_BASIC_INTERFACE_H__
