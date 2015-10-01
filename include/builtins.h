#ifndef __DASHMM_BUILT_INS_H__
#define __DASHMM_BUILT_INS_H__


#include "include/expansion.h"
#include "include/method.h"


namespace dashmm {


Expansion *laplace_cartesian_expansion();
Expansion *laplace_spherical_harmonic_expansion(int n_digits);

Method *bh_method(double theta);
Method *fmm_classic_method();


//TODO: add a few routines that will provide the MethodDesc or ExpansionDesc
// for these builtins. This will allow the user to make modified versions of
// these easily. Also, we can use these routines during init to register the
// built-ins ourselves. This means the implementation of these manually
// constructs the description each time it is called, but whatever. That is
// init code.


} // namespace dashmm


#endif // __DASHMM_BUILT_INS_H__
