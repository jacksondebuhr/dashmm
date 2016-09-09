#ifndef __DASHMM_BROADCAST_H__
#define __DASHMM_BROADCAST_H__


#include "hpx/hpx.h"


namespace dashmm {


extern hpx_action_t broadcast_value_action;


/// This broadcasts the given value from rank 0 to the other ranks
template <typename T>
void broadcast(T *value) {
  hpx_run(&broadcast_value_action, value, value, sizeof(T));
}


} // dashmm


#endif // __DASHMM_BROADCAST_H__