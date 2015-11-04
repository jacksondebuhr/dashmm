#include "include/builtins.h"

// C/C++

#include <hpx/hpx.h>

#include "include/bh_method.h"
#include "include/expansion.h"
#include "include/laplace_com.h"
#include "include/method.h"



namespace dashmm {


Method *bh_method(double theta) {
  Method *retval = new BHMethod{theta};
  return retval;
}


Method *bh_method_deserialize(size_t size, MethodSerial *data) {
  assert(size == (sizeof(MethodSerial) + sizeof(double)));
  assert(data->size == sizeof(double));
  double *theta = static_cast<double *>(data->data);
  return bh_method(*theta);
}
HPX_ACTION(HPX_FUNCTION, 0,
           bh_method_deserialize_action, bh_method_deserialize,
           HPX_SIZE_T, HPX_POINTER);


void register_built_in_methods() {
  assert(register_method(kMethodBH, bh_method_deserialize_action));
  //more here...
}


Expansion *laplace_COM_expansion() {
  Expansion *retval = new LaplaceCOM{Point{0.0, 0.0, 0.0}};
  return retval;
}


Expansion *laplace_COM_deserialize(size_t size, ExpansionSerial *data) {
  assert(size == (sizeof(ExpansionSerial) + sizeof(double) * 10));
  assert(data->size == 10 * sizeof(double));
  LaplaceCOM *retval = new LaplaceCOM{Point{0.0, 0.0, 0.0}};
  double *vals = static_cast<double *>(data->data);
  retval->set_mtot(vals[0]);
  retval->set_xcom(&vals[1]);
  retval->set_Q(&vals[4]);
  return retval;
}
HPX_ACTION(HPX_FUNCTION, 0,
           laplace_COM_deserialize_action, laplace_COM_deserialize,
           HPX_SIZE_T, HPX_POINTER);


void register_built_in_expansions() {
  assert(register_expansion(kExpansionLaplaceCOM,
                            laplace_COM_deserialize_action));
  //more here
}


}  // namespace dashmm
