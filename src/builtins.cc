/// \file src/builtins.cc
/// \brief Implementation of built-in sources and expansions


#include "include/builtins.h"

#include <hpx/hpx.h>

#include "include/bh_method.h"
#include "include/direct_method.h"
#include "include/expansion.h"
#include "include/laplace_com.h"
#include "include/method.h"



namespace dashmm {


Method *bh_method(double theta) {
  Method *retval = new BHMethod{theta};
  return retval;
}


Method *bh_method_create_handler(size_t size, MethodSerial *data) {
  assert(size == (sizeof(MethodSerial) + sizeof(double)));
  assert(data->size == sizeof(double));
  double *theta = static_cast<double *>(data->data);
  return bh_method(*theta);
}
HPX_ACTION(HPX_FUNCTION, 0,
           bh_method_create_action, bh_method_create_handler,
           HPX_SIZE_T, HPX_POINTER);


Method *direct_method() {
  Method *retval = new DirectMethod{};
  return retval;
}


Method *direct_method_create_handler(size_t size, MethodSerial *data) {
  assert(size == sizeof(MethodSerial));
  return direct_method();
}
HPX_ACTION(HPX_FUNCTION, 0,
           direct_method_create_action, direct_method_create_handler,
           HPX_SIZE_T, HPX_POINTER);


void register_built_in_methods() {
  assert(kSuccess == register_method(kMethodBH, bh_method_create_action, 0));
  assert(kSuccess == register_method(
                        kMethodDirect, direct_method_create_action, 0));
}


Expansion *laplace_COM_expansion() {
  Expansion *retval = new LaplaceCOM{Point{0.0, 0.0, 0.0}};
  return retval;
}


Expansion *laplace_COM_create_handler(double x, double y, double z) {
  LaplaceCOM *retval = new LaplaceCOM{Point{x, y, z}};
  return retval;
}
HPX_ACTION(HPX_FUNCTION, 0,
           laplace_COM_create_action, laplace_COM_create_handler,
           HPX_DOUBLE, HPX_DOUBLE, HPX_DOUBLE);


Expansion *laplace_COM_interpret_handler(void *data, size_t size) {
  LaplaceCOM *retval = new LaplaceCOM{(LaplaceCOMData *)data, size};
  return retval;
}
HPX_ACTION(HPX_FUNCTION, 0,
           laplace_COM_interpret_action, laplace_COM_interpret_handler,
           HPX_SIZE_T, HPX_POINTER);


void register_built_in_expansions() {
  assert(kSuccess == register_expansion(kExpansionLaplaceCOM,
                                        laplace_COM_create_action,
                                        laplace_COM_interpret_action, 0));
  //more here
}


}  // namespace dashmm
