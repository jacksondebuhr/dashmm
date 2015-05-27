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
/// This is the handle to the built-in Yukawa kernel.
///
extern dashmm_handle_t DASHMM_YUKAWA_POTENTIAL;


// ==================================================================
// Basic Interface to DASHMM
// ==================================================================


///
/// \brief Initializes the runtime and DASHMM
///
/// JD: We will need to speak with the developers about how we will need to 
///     initialize the runtime and so on, for various use cases
int dashmm_init(...)


///
/// \brief Cleans up DASHMM and halts the runtime
///
/// JD: Again, we shall need to speak with the devs to see what this routine
///     may need.
int dashmm_finalize(...)


///
/// \brief Evaluate forces and/or potentials using a 
///        hierarchical multipole method
///
/// This function evalutates the potential and/or the force for a given 
/// potential and using a given evaluation method. The selected method and
/// kernel can be either a built-in method or kernel or one specified by the
/// user. DASHMM supports evaluation in either a FMM-like or BH-like manner,
/// so any user specified methods will be variations on one of these two
/// overall classes of methods.
///
/// Depending on the class of method employed, @param accuracy will set the 
/// number of digits of accuracy (FMM-like) or the number of terms in the 
/// multipole expansion (BH-like).
///
/// The source and target points are specified by providing the handle to a 
/// globally allocated buffer that is registered with DASHMM. In addition to 
/// the handle to the buffer, this routine needs the offset in the records of
/// the relevant buffer of the positions. For the source points, the user
/// must also specify the offset in the records to the charges of the source
/// points. For the target locations, the user must specify the offset to the
/// results in the record. 
///
/// Users can use the same global buffer for the source and target points. 
/// This routine will not modify the input buffer, with the single exception
/// that if the source and target buffers are the same, the results will be 
/// stored in the correct location.
///
/// This routine will compute either the potential at the target locations,
/// the force at the target location, or both. This behavior is selected
/// with the @param result_selection parameter.
///
/// Finally, the user may optionally provide the address of a dashmm_stats_t
/// object in which to store collected timing information.
///
/// \param method - the method that will be used for the evaluation; this can
///                 be one of the built in methods, or can be defined by the
///                 user.
/// \param kernel - the kernel that will be used for the evaluation; this can
///                 be one of the built in kernels, or can be a kernel defined
///                 by the user.
/// \param accuracy - the parameter controls the accuracy of the chosen method;
///                 giving either the requested accuracy for 
///                 FMM-like methods, or the opening angle for BH-like methods.
/// \param source - the handle to the global buffer holding the source data.
///                 this routine will not modify the input data, unless the
///                 same global buffers are used as @param source and 
///                 @param target.
/// \param source_position_offset - offset in each record of @param source to
///                 the source positions
/// \param source_charge_offset - offset in each record @param source to the
///                 charges
/// \param result_selection - an ORed together list of DASHMM_POTENTIAL and
///                 DASHMM_FORCE so select if potentials, forces or both are
///                 to be computed. If none are selected, this function will
///                 perform no work, and return DASHMM_DOMAIN_ERROR.
/// \param target - the handle to the global buffer holding the target data
/// \param target_position_offset - the offset in each record of @param target
///                 to the position of the target points.
/// \param target_result_offset - the offset in each record of @param target to
///                 the place to store the results of the evaluation. There will
///                 be between 1 and 4 values stored, depending on
///                 @param result_selection, with the potential stored first,
///                 followed by the three components of the force.
/// \param stats - the address of a dashmm_stats_t to collect timing statistics
///                 for the evaluation, or NULL if timing information is not 
///                 needed.
///
/// \return - DASHMM_SUCCESS on successful evaluation; or other values on error.
/// 
int dashmm_evaluate(dashmm_handle_t method,
                    dashmm_handle_t kernel,
                    double accuracy,
                    dashmm_handle_t source, 
                    size_t source_position_offset,
                    size_t source_charge_offset,
                    int result_selection,
                    dashmm_handle_t target,
                    size_t target_position_offset,
                    size_t target_result_offset,
                    dashmm_stats_t *stats);


///
/// \brief Allocate an array in the global address space
///
/// This routine allocated a buffer in the global address space provided by
/// DASHMM. The buffer is arranged into a series of records with a fixed size.
/// These records are distributed among the localities available to the 
/// DASHMM system. For more control over the layout of these records please 
/// see <some other function that is not defined yet>.
///
/// The returned handle can be used in all DASHMM routines to refer to the
/// allocated array. In particular, the source and target data must reside
/// in the global address space before dashmm_evaluate may be used.
///
/// \param records     - the number of records to allocate
/// \param record_size - the size in bytes of each record
/// 
/// \return - the handle to the allocated memory; if this value is 
//            INVALID_HANDLE, there was an error during the allocation
///
dashmm_handle_t dashmm_array_alloc(uint64_t records, size_t record_size);


///
/// \brief Free an object in the global address space
///
/// The releases global memory previously allocated by DASHMM.
/// @param handle must be a handle to memory that has been allocated and
/// registered with DASHMM, otherwise this routine will return an error.
/// Additionally, it is an error to attempt to free system-defined objects.
///
/// \param handle - the handle to the object in question
///
/// \return - DASHMM_DOMAIN_ERROR if handle is not a user array object;
///           DASHMM_RUNTIME_ERROR if an error from the runtime is encountered;
///           DASHMM_SUCCESS otherwise
///
int dashmm_array_free(dashmm_handle_t handle);


///
/// \brief Write memory into a global buffer managed by DASHMM
///
/// This routine writes local memory into a global buffer registerd with the
/// DASHMM system. This writes @param length bytes starting at @param offset
/// bytes from the start of @param record in the global buffer identified by
/// @param handle. @param data should point to the local memory that is to 
/// be written into the global buffer.
///
/// \param handle - the handle to the global object
/// \param record - selects which record to address
/// \param offset - the offset in bytes from the start of the record to begin
///                 writing
/// \param length - the number of bytes to write
/// \param data   - the address of the data to write into the global buffer
///
/// \return - DASHMM_SUCCESS or DASHMM_DOMAIN_ERROR if record is larger than
///           the number of records in the provided array.
///
int dashmm_array_memput(dashmm_handle_t handle, 
                        uint64_t record, 
                        size_t offset,
                        size_t length,
                        void *data);


///
/// \brief Read memory from a global buffer managed by DASHMM
///
/// This routine reads from a global buffer registerd with the DASHMM system
/// into local memory. This reads @param length bytes starting at @param offset
/// bytes from the start of @param record in the global buffer identified by
/// @param handle. @param data should point to the local memory that will hold 
/// the data read from global memory.
///
/// \param handle - the handle to the global object
/// \param record - selects which record to address
/// \param offset - the offset in bytes from the start of the record to begin
///                 reading
/// \param length - the number of bytes to read
/// \param data   - address of local memory to hold the data read from global 
///                 memory
///
/// \return - DASHMM_SUCCESS or DASHMM_DOMAIN_ERROR if record is larger than
///           the number of records in the provided array.
///
int dashmm_array_memget(dashmm_handle_t handle,
                        uint64_t record,
                        size_t offset,
                        size_t length,
                        void *data);


// ==================================================================
// Advanced Interface to DASHMM
// ==================================================================
//      The intent is that the user will need only to include a single file
//      to use the library. As such, the full interface will appear in this
//      header.
//      Of course, this could just mean that this includes a bunch of other
//      files.


#endif // __DASHMM_INTERFACE_H__
