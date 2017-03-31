// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2017, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license. See the LICENSE file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================


#ifndef __DASHMM_TRACE_EVENTS_H__
#define __DASHMM_TRACE_EVENTS_H__


/// \file
/// \brief Definitions for DASHMM event traces


#ifdef DASHMM_INSTRUMENTATION

#include <hpx/builtins.h>

#define ENABLE_INSTRUMENTATION
#include <libhpx/instrumentation.h>


typedef enum {
#define LIBHPX_EVENT(class, event, ...) TRACE_EVENT_##class##_##event,
#define LIBHPX_REMOVE_EVENT(class, event, ...) TRACE_EVENT_##class##_##event,
# include <libhpx/events.def>
#undef LIBHPX_REMOVE_EVENT
#undef LIBHPX_EVENT
} libhpx_trace_events_t;


# define _DECL0() void
# define _DECL2(t0,n0) t0 u0
# define _DECL4(t0,n0,t1,n1) t0 u0, t1 u1
# define _DECL6(t0,n0,t1,n1,t2,n2) t0 u0, t1 u1, t2 u2
# define _DECL8(t0,n0,t1,n1,t2,n2,t3,n3) t0 u0, t1 u1, t2 u2, t3 u3
# define _DECL10(t0,n0,t1,n1,t2,n2,t3,n3,t4,n4) t0 u0, t1 u1, t2 u2, t3 u3, t4 u4
# define _DECLN(...) _HPX_CAT2(_DECL, __HPX_NARGS(__VA_ARGS__))(__VA_ARGS__)
# define _ARGS0()
# define _ARGS2(t0,n0) , u0
# define _ARGS4(t0,n0,t1,n1) , u0, u1
# define _ARGS6(t0,n0,t1,n1,t2,n2) , u0, u1, u2
# define _ARGS8(t0,n0,t1,n1,t2,n2,t3,n3) , u0, u1, u2, u3
# define _ARGS10(t0,n0,t1,n1,t2,n2,t3,n3,t4,n4) , u0, u1, u2, u3, u4
# define _ARGSN(...) _HPX_CAT2(_ARGS, __HPX_NARGS(__VA_ARGS__))(__VA_ARGS__)
# define LIBHPX_EVENT(class, event, ...)                                \
  static inline void                                                    \
  EVENT_##class##_##event(_DECLN(__VA_ARGS__)) {                        \
    trace_append(HPX_TRACE_##class, TRACE_EVENT_##class##_##event       \
                 _ARGSN(__VA_ARGS__));                                  \
  }
# define LIBHPX_REMOVE_EVENT(class, event, ...)                         \
  static inline void                                                    \
  EVENT_##class##_##event(_DECLN(__VA_ARGS__)) {}
# include <libhpx/events.def>
# undef LIBHPX_REMOVE_EVENT
# undef LIBHPX_EVENT
# undef _ARGS0
# undef _ARGS2
# undef _ARGS4
# undef _ARGS6
# undef _ARGS8
# undef _ARGS10
# undef _ARGSN
# undef _DECL0
# undef _DECL2
# undef _DECL4
# undef _DECL6
# undef _DECL8
# undef _DECL10
# undef _DECLN

#else

#define EVENT_TRACE_DASHMM_STOT_BEGIN(...)
#define EVENT_TRACE_DASHMM_STOT_END(...)
#define EVENT_TRACE_DASHMM_MTOT_BEGIN(...)
#define EVENT_TRACE_DASHMM_MTOT_END(...)
#define EVENT_TRACE_DASHMM_LTOT_BEGIN(...)
#define EVENT_TRACE_DASHMM_LTOT_END(...)
#define EVENT_TRACE_DASHMM_STOM_BEGIN(...)
#define EVENT_TRACE_DASHMM_STOM_END(...)
#define EVENT_TRACE_DASHMM_STOL_BEGIN(...)
#define EVENT_TRACE_DASHMM_STOL_END(...)
#define EVENT_TRACE_DASHMM_ELCO_BEGIN(...)
#define EVENT_TRACE_DASHMM_ELCO_END(...)
#define EVENT_TRACE_DASHMM_MTOM_BEGIN(...)
#define EVENT_TRACE_DASHMM_MTOM_END(...)
#define EVENT_TRACE_DASHMM_MTOL_BEGIN(...)
#define EVENT_TRACE_DASHMM_MTOL_END(...)
#define EVENT_TRACE_DASHMM_LTOL_BEGIN(...)
#define EVENT_TRACE_DASHMM_LTOL_END(...)
#define EVENT_TRACE_DASHMM_MTOI_BEGIN(...)
#define EVENT_TRACE_DASHMM_MTOI_END(...)
#define EVENT_TRACE_DASHMM_ITOI_BEGIN(...)
#define EVENT_TRACE_DASHMM_ITOI_END(...)
#define EVENT_TRACE_DASHMM_ITOL_BEGIN(...)
#define EVENT_TRACE_DASHMM_ITOL_END(...)
#define EVENT_TRACE_DASHMM_ZEROREF(...)

#endif


#endif