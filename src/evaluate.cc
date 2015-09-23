//C++ stuff

#include "include/expansion.h"
#include "include/method.h"
#include "include/node.h"
#include "include/types.h"


namespace dashmm {


/////////////////////////////////////////////////////////////////////
// Actions
/////////////////////////////////////////////////////////////////////


int evaluate_handler(hpx_addr_t sources, hpx_addr_t targets,
                     int refinement_limit,
                   /*TODO something here for method and expansion*/) {
  //find domain bounds
}
HPX_ACTION();


/////////////////////////////////////////////////////////////////////
// Interface
/////////////////////////////////////////////////////////////////////


int evaluate(ObjectHandle sources, ObjectHandle targets, int refinement_limit,
             Method *method, Expansion *expansion) {
  if (!method->compatible_with(expansion)) {
    return kIncompatible;
  }

  //find domain bounds

  //build trees / do the work

  //delete trees

  return kSuccess;
}


} // namespace dashmm
