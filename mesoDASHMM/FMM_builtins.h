#ifndef __FMM_BUILTINS_H__
#define __FMM_BUILTINS_H__


#include "expansion.h"
#include "method.h"


namespace dashmm {


Method *fmm_method();
Expansion *fmm_expansion(const int p);

//Again, these will eventially be hidden away from the user.
void fmm_builtins_init();
void fmm_builtins_fini();


} //namespace dashmm


#endif
