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

#ifndef __DASHMM_BASIC_TYPES_H__
#define __DASHMM_BASIC_TYPES_H__


#include "hpx/hpx.h"


///
/// \type dashmm_status_t
///
/// Enumeration of the possible return codes from DASHMM library calls.
///
typedef enum {
  DASHMM_ERROR = -1;
  DASHMM_SUCCESS = 0;
  DASHMM_DOMAIN_ERROR = 1;
  //and so on...
} dashmm_status_t;


///
/// \type dashmm_result_options_t
///
/// This enumeration is used to specify options effecting the results produced
/// byt DASHMM evaluations.
///
typedef enum {
  DASHMM_POTENTIAL = 1;
  DASHMM_FORCE     = 2;
  //etc?...
} dashmm_result_options_t;


///
/// \type dashmm_handle_t
///
/// The user-presented handle to DASHMM objects. As can be seen, this is just
/// a renamed GAS address from HPX-5.
typedef hpx_addr_t dashmm_handle_t;

#define INVALID_HANDLE HPX_NULL


///
/// \type dashmm_stats_t
///
/// This type holds basic timing statistics for an evaluation using DASHMM,
/// including the total time take, the min and max times over each locality in
/// the system, the median times and the variance.
///
typedef struct {
  double t_total;
  double t_min;
  double t_max;
  double t_mean;
  double t_variance;
} dashmm_stats_t;


#endif
