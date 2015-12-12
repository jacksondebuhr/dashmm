//C++ stuff

#include "hpx/hpx.h"

#include "include/builtins.h"
#include "include/types.h"


namespace dashmm {


/////////////////////////////////////////////////////////////////////
// Actions
/////////////////////////////////////////////////////////////////////


int init_handler(void *UNUSED, size_t UNWANTED) {
  //do a broadcast of some stuff, who knows...

  //Create the registration tables
  init_method_table();
  init_expansion_table();

  //These register any built-in methods or expansions.
  register_built_in_methods();
  register_built_in_expansions();

  hpx_exit(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED, init_action, init_handler,
           HPX_POINTER, HPX_SIZE_T);


int fini_handler(void *UNUSED, size_t UNWANTED) {
  //do a broadcast of something?

  //Clear out the registrations
  fini_method_table();
  fini_expansion_table();

  hpx_exit(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED, fini_action, fini_handler,
           HPX_POINTER, HPX_SIZE_T);


/////////////////////////////////////////////////////////////////////
// Interface
/////////////////////////////////////////////////////////////////////


ReturnCode init(int *argc, char ***argv) {
  if (HPX_SUCCESS != hpx_init(argc, argv)) {
    return kRuntimeError;
  }

  //now run the initialization action
  if (HPX_SUCCESS != hpx_run(&init_action, nullptr, 0)) {
    return kInitError;
  }

  return kSuccess;
}


ReturnCode finalize() {
  //run the finalization action
  if (HPX_SUCCESS != hpx_run(&fini_action, nullptr, 0)) {
    return kFiniError;
  }

  //shutdown the runtime
  hpx_finalize();

  return kSuccess;
}


} // namespace dashmm
