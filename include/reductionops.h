#ifndef __DASHMM_REDUCTION_OPS_H__
#define __DASHMM_REDUCTION_OPS_H__


/// \file include/reductionops.h
/// \brief Action identifiers for common reduction operations


#include <hpx/hpx.h>


namespace dashmm {


/// Identity operation for integer summation
extern hpx_action_t int_sum_ident_op;

/// Operation for integer summation
extern hpx_action_t int_sum_op;


} // namespace dashmm


#endif // __DASHMM_REDUCTION_OPS_H__
