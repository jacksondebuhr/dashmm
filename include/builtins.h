#ifndef __DASHMM_BUILT_INS_H__
#define __DASHMM_BUILT_INS_H__


#include "include/expansion.h"
#include "include/method.h"


namespace dashmm {


Expansion *laplace_cartesian_expansion();

Expansion *laplace_spherical_harmonic_expansion(int n_digits);

Method *bh_method(double theta);

Method *fmm_classic_method();


} // namespace dashmm


#endif // __DASHMM_BUILT_INS_H__
