#ifndef __DASHMM_BUILT_INS_H__
#define __DASHMM_BUILT_INS_H__


#include "include/expansion.h"
#include "include/method.h"


namespace dashmm {


//These provide instances of these expansions
Expansion *laplace_COM_expansion();
Expansion *laplace_sph_expansion(int n_digits);
Expansion *laplace_sph_exp_expansion(int n_digits);

//These provide instances of these methods
Method *bh_method(double theta);
Method *fmm_method();
Method *fmm_exp_method();


} // namespace dashmm


#endif // __DASHMM_BUILT_INS_H__
