//C++ stuff

#include "hpx/hpx.h"

#include "include/types.h"


namespace dashmm {


/////////////////////////////////////////////////////////////////////
// Actions
/////////////////////////////////////////////////////////////////////


int init_handler(void *UNUSED, size_t UNWANTED) {
  //do a broadcast of some stuff, who knows...
  hpx_exit(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED, init_action, init_handler,
           HPX_POINTER, HPX_SIZE_T);


int fini_handler(void *UNUSED, size_t UNWANTED) {
  //do a broadcast of something?
  hpx_exit(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED, fini_action, fini_handler,
           HPX_POINTER, HPX_SIZE_T);


/////////////////////////////////////////////////////////////////////
// Interface
/////////////////////////////////////////////////////////////////////


int init(int *argc, char ***argv) {
  if (HPX_SUCCESS != hpx_init(argc, argv)) {
    return kRuntimeError;
  }

  //now run the initialization action
  if (HPX_SUCCESS != hpx_run(init_action, nullptr, 0)) {
    return kInitError;
  }

  return kSuccess;
}


int finalize() {
  //run the finalization action
  if (HPX_SUCCESS != hpx_run(fini_action, nullptr, 0)) {
    return kFiniError;
  }

  //shutdown the runtime
  hpx_finalize();

  return kSuccess;
}


} // namespace dashmm
