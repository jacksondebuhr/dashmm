#ifndef __DASHMM_REDUCTION_OPS_H__
#define __DASHMM_REDUCTION_OPS_H__


#include <hpx/hpx.h>


namespace dashmm {


//Reduction operators for a sum of n integers
extern hpx_action_t int_sum_ident_op;
extern hpx_action_t int_sum_op;


} // namespace dashmm


#endif // __DASHMM_REDUCTION_OPS_H__
