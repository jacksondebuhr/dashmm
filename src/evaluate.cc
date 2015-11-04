//C++ stuff

#include "include/expansion.h"
#include "include/method.h"
#include "include/node.h"
#include "include/types.h"


namespace dashmm {


/////////////////////////////////////////////////////////////////////
// Utility routines
/////////////////////////////////////////////////////////////////////


hpx_addr_t pack_sources(hpx_addr_t user_data, int pos_offset, int q_offset,
                        int *nparts) {
  //NOTE: SMP assumption
  ArrayMetaData *meta{nullptr};
  assert(hpx_gas_try_pin(user_data, (void **)&meta));
  *nparts = meta->count;

  char *user{nullptr};
  assert(hpx_gas_try_pin(meta->data, (void **)&user));

  hpx_addr_t retval = hpx_gas_alloc_local_at_sync(1, meta->count * meta->size,
                                                  0, HPX_HERE);
  if (retval != HPX_NULL) {
    Source *sources{nullptr};
    assert(hpx_gas_try_pin(retval, (void **)&sources));

    //TODO some kind of placement new here...
    for (size_t i = 0; i < meta->count; ++i) {
      //
    }
  }

  return retval;
}


hpx_addr_t pack_targets(hpx_addr_t user_data, int pos_offset, int *nparts) {
  //
}


void unpack_targets(hpx_addr_t results, hpx_addr_t user_data, int phi_offset) {
  //
}


/////////////////////////////////////////////////////////////////////
// Actions
/////////////////////////////////////////////////////////////////////


struct EvaluateParams {
  hpx_addr_t sources;
  int spos_offset;
  int q_offset;
  hpx_addr_t targets;
  int tpos_offset;
  int phi_offset;
  int refinement_limit;
  size_t method_size;
  size_t expansion_size;
  char data[];
};

//TODO: Go back through to see if there is some overlap we can perform on
// these things at all.
int evaluate_handler(EvaluateParams *parms, size_t total_size) {
  //create our method and expansion from the parameters
  MethodSerial *method_serial = static_cast<MethodSerial *>(parms->data);
  std::unique_ptr<Method> method =
      create_method(method_serial->type, parms->method_size, method_serial);

  ExpansionSerial *expansion_serial =
      static_cast<ExpansionSerial *>(parms->data + parms->method_size);
  std::unique_ptr<Expansion> expansion = create_expansion(
        expansion_serial->type, parms->expansion_size, expansion_serial);

  //copy user data into internal data
  int n_sources{0};
  hpx_addr_t source_parts = pack_sources(parms->sources, parms->spos_offset,
                                         parms->q_offset, &n_sources);
  int n_targets{0};
  hpx_addr_t target_parts = pack_targets(parms->targets, parms->tpos_offset,
                                         &n_targets);

  //find domain bounds

  //build trees/do work

  //delete trees

  //copy results back into user data
  unpack_targets(target_parts, parms->targets, parms->phi_offset);

  //free intermediate data
  hpx_gas_free_sync(source_parts);
  hpx_gas_free_sync(target_parts);

  //return
  return HPX_SUCCESS;
}
HPX_ACTION(HPX_ACTION, HPX_MARSHALLED,
           evaluate_action, evaluate_handler,
           HPX_POINTER, HPX_SIZE_T);


/////////////////////////////////////////////////////////////////////
// Interface
/////////////////////////////////////////////////////////////////////


ReturnCode evaluate(ObjectHandle sources, int spos_offset, int q_offset,
                    ObjectHandle targets, int tpos_offset, int phi_offset,
                    int refinement_limit,
                    Method *method, Expansion *expansion) {
  if (!method->compatible_with(expansion)) {
    return kIncompatible;
  }

  //Marshall our arguments
  MethodSerialPtr method_serial = method->serialize(false);
  ExpansionSerialPtr expansion_serial = expansion->serialize(false);
  size_t method_size = sizeof(MethodSerial) + method_serial->size;
  size_t expansion_size = sizeof(ExpansionSerial) + expansion_serial->size;
  size_t total_size = method_size + expansion_size + sizeof(EvaluateParms);
  EvaluateParms *args = static_cast<EvaluateParms *>(malloc(total_size));
  assert(args);
  args->sources = sources;
  args->spos_offset = spos_offset;
  args->q_offset = q_offset;
  args->targets = targets;
  args->tpos_offset = tpos_offset;
  args->phi_offset = phi_offset;
  args->refinement_limit = refinement_limit;
  args->method_size = method_size;
  args->expansion_size = expansion_size;
  memcpy(data, method_serial.get(), method_size);
  memcpy(data + method_size, expansion_serial.get(), expansion_size);

  if (HPX_SUCCESS != hpx_run(evaluate_action, args, total_size)) {
    return kRuntimeError;
  }

  return kSuccess;
}


} // namespace dashmm
