#ifndef __DASHMM_BUILT_INS_H__
#define __DASHMM_BUILT_INS_H__


#include "include/types.h"


namespace dashmm {


//expansion for laplace using cartesian coordianates, and which includes up to
// quadrupole terms, and which is expanded about the center of mass of the box.
// This really should not be used with charges of both sign.
//TODO: decide if we should extend to have one that works with both signs?
// This would not be too hard. It would just use the center of the box as the
// expansion location instead of the COM. Unless we get lucky, there is likely
// no place in the both signs case that will give a zero dipole moment.
//So perhaps we just go ahead and use the center?
//Or should we just make those things arguments to this call, so expansion
// location (center or center of mass). Using com would requrie
Expansion *laplace_com_expansion();

Expansion *laplace_sph_expansion(int n_digits);

Method *bh_method(double theta);

Method *fmm_classic_method();


} // namespace dashmm


#endif
