#ifndef __DASHMM_BUILT_INS_H__
#define __DASHMM_BUILT_INS_H__


/// \file include/builtins.h
/// \brief Builtin methods and expansions


#include "include/expansion.h"
#include "include/method.h"


namespace dashmm {


/// Provide a Laplace Center of Mass expansion
///
/// This creates a new object, the caller assumes ownership of the returned
/// object.
///
/// \returns - the expansion; caller assumes ownership
Expansion *laplace_COM_expansion();

/// Provide a Laplace spherical harmonic expansion
///
/// This creates a new object, the caller assumes ownership of the returned
/// object.
///
/// \returns - the expansion; caller assumes ownership
Expansion *laplace_sph_expansion(int n_digits);

/// Provide a Laplace spherical harmonic with exponential operators expansion
///
/// This creates a new object, the caller assumes ownership of the returned
/// object.
///
/// \returns - the expansion; caller assumes ownership
Expansion *laplace_sph_exp_expansion(int n_digits);

/// Provide a BH method object
///
/// This creates a new object, the caller assumes ownership of the returned
/// object. This method is the classic Barnes-Hut method.
///
/// \returns - the method; caller assumes ownership
Method *bh_method(double theta);

/// Provide an FMM method object
///
/// This creates a new object, the caller assumes ownership of the returned
/// object. This method is the classic Fast Multipole Method.
///
/// \returns - the method; caller assumes ownership
Method *fmm_method();

/// Provide an FMM method object that uses exponential operators
///
/// This creates a new object, the caller assumes ownership of the returned
/// object. This method is the Fast Multipole Method that uses Merge-and-shift
/// to reduce the arithmetic complexity.
///
/// \returns - the method; caller assumes ownership
Method *fmm_exp_method();


// NOTE: Not intended for users

/// Register the built in methods with DASHMM
///
/// This method is not intended for end-users.
void register_built_in_methods();

/// Register the built in expansions with DASHMM
///
/// This method is not intended for end-users.
void register_built_in_expansions();


} // namespace dashmm


#endif // __DASHMM_BUILT_INS_H__
