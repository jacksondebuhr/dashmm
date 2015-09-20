#ifndef __BH_BUILTINS_H__
#define __BH_BUILTINS_H__


#include "expansion.h"
#include "method.h"


namespace dashmm {


//These are part of the interface
Method *bh_method(double crit);
Expansion *bh_expansion();

//These are not, but are provided for convenience, and so we do not have to
// expose the globals
//In the future, we would likely put these declarations in a separate header
// and would only include that in the dashmm.cc file .
void bh_builtins_init();
void bh_builtins_fini();


}


#endif
