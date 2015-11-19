//C++ stuff

#include "include/expansion.h"
#include "include/method.h"
#include "include/node.h"
#include "include/types.h"


namespace dashmm {


/////////////////////////////////////////////////////////////////////
// Utility routines
/////////////////////////////////////////////////////////////////////


DomainGeometry cubify_domain(hpx_addr_t source_bounds,
                             hpx_addr_t target_bounds) {
  double s_bounds[6];
  hpx_lco_get(source_bounds, sizeof(double) * 6, s_bounds);

  double t_bounds[6];
  hpx_lco_get(target_bounds, sizeof(double) * 6, t_bounds);

  //TODO Can't I use the DomainGeometry Constructor here?
  //cubify the domain
  Point low{s_bounds[0] < t_bounds[0] ? s_bounds[0] : t_bounds[0],
            s_bounds[1] < t_bounds[1] ? s_bounds[1] : t_bounds[1],
            s_bounds[2] < t_bounds[2] ? s_bounds[2] : t_bounds[2]};
  Point high{s_bounds[3] > t_bounds[3] ? s_bounds[3] : t_bounds[3],
             s_bounds[4] > t_bounds[4] ? s_bounds[4] : t_bounds[4],
             s_bounds[5] > t_bounds[5] ? s_bounds[5] : t_bounds[5]};
  Point center = point_add(low.scale(0.5), high.scale(0.5));
  Point sizes = point_sub(high, low);
  double size = sizes.x() > sizes.y() ? sizes.x() : sizes.y();
  size = size > sizes.z() ? size : sizes.z();
  size *= 0.5001;
  Point offset{-size, -size, -size};
  low = point_add(high, offset);
  size *= 2.0;

  return DomainGeometry{low, size};
}


/////////////////////////////////////////////////////////////////////
// Actions
/////////////////////////////////////////////////////////////////////


struct PackDataResult {
  hpx_addr_t packed;
  int count;
};


//TODO make this continue a sourceref
int pack_sources_handler(hpx_addr_t user_data, int pos_offset, int q_offset) {
  //NOTE: SMP assumptions all over this function
  ArrayMetaData *meta{nullptr};
  assert(hpx_gas_try_pin(user_data, (void **)&meta));

  char *user{nullptr};
  assert(hpx_gas_try_pin(meta->data, (void **)&user));

  PackDataResult retval{};
  retval.packed =
      hpx_gas_alloc_local_at_sync(1, meta->count * sizeof(Source), 0, HPX_HERE);
  retval.count = meta->count;

  if (retval.packed != HPX_NULL) {
    Source *sources{nullptr};
    assert(hpx_gas_try_pin(retval.packed, (void **)&sources));

    for (size_t i = 0; i < meta->count; ++i) {
      char *record_base = user[i * meta->size];
      double *pos = static_cast<double *>(record_base + pos_offset);
      double *q = static_cast<double *>(record_base + q_offset);
      sources[i]->set_position(Point{pos[0], pos[1], pos[2]});
      sources[i]->set_charge(*q);
    }

    hpx_gas_unpin(retval.packed);
  }

  hpx_gas_unpin(meta->data);
  hpx_gas_unpin(user_data);

  HPX_THREAD_CONTINUE(retval);
}
HPX_ACTION(HPX_DEFAULT, 0,
           pack_sources_action, pack_sources_handler,
           HPX_ADDR, HPX_INT, HPX_INT);


int pack_targets_handler(hpx_addr_t user_data, int pos_offset) {
  //NOTE: SMP assumptions
  ArrayMetaData *meta{nullptr};
  assert(hpx_gas_try_pin(user_data, (void **)&meta));

  char *user{nullptr};
  assert(hpx_gas_try_pin(meta->data, (void **)&user));

  PackDataResult retval{};
  retval.packed =
      hpx_gas_alloc_local_at_sync(1, meta->count * sizeof(Target), 0, HPX_HERE);
  retval.count = meta->count;
  if (retval.packed != HPX_NULL) {
    Target *targets{nullptr};
    assert(hpx_gas_try_pin(retval, (void **)&targets));

    for (size_t i = 0; i < meta->count; ++i) {
      char *record_base = user[i * meta->size];
      double *pos = static_cast<double *>(record_base + pos_offset);
      targets[i]->set_position(Point{pos[0], pos[1], pos[2]});
      targets[i]->set_index(i);
    }

    hpx_gas_unpin(retval);
  }

  hpx_gas_unpin(meta->data);
  hpx_gas_unpin(user_data);

  HPX_THREAD_CONTINUE(retval);
}
HPX_ACTION(HPX_DEFAULT, 0,
           pack_targets_action, pack_targets_handler,
           HPX_ADDR, HPX_INT);


int unpack_targets_handler(hpx_addr_t results, hpx_addr_t user_data, int phi_offset) {
  //NOTE: SMP assumptions
  ArrayMetaData *meta{nullptr};
  assert(hpx_gas_try_pin(user_data, (void **)meta));

  char *user{nullptr};
  assert(hpx_gas_try_pin(meta->data, (void **)&user));

  Target *targets{nullptr};
  assert(hpx_gas_try_pin(results, (void **)&targets));

  for (size_t i = 0; i < meta->count; ++i) {
    size_t idx = targets[i].index();
    char *record_base = user[idx * meta->size];
    double *phi = static_cast<double *>(record_base + phi_offset);
    *phi = targets[i].phi();
  }

  hpx_gas_unpin(results);
  hpx_gas_unpin(meta->data);
  hpx_gas_unpin(user_data);

  return HPX_SUCCESS;
}
HPX_ACTION(HPX_DEFAULT, 0,
           unpack_targets_action, unpack_targets_handler,
           HPX_ADDR, HPX_ADDR, HPX_INT);


int find_source_domain_handler(Source *sources, int n_sources) {
  //NOTE: SMP assumptions

  //These are the three low bounds, followed by the three high bounds
  double bounds[6]{1.0e34, 1.0e34, 1.0e34, -1.0e34, -1.0e34, -1.0e34};

  for (int i = 0; i < n_sources; ++i) {
    Point p = sources[i].position();
    if (p.x() < bounds[0]) bounds[0] = p.x();
    if (p.x() > bounds[3]) bounds[3] = p.x();
    if (p.y() < bounds[1]) bounds[1] = p.y();
    if (p.y() > bounds[4]) bounds[4] = p.y();
    if (p.z() < bounds[2]) bounds[2] = p.z();
    if (p.z() > bounds[5]) bounds[5] = p.z();
  }

  hpx_thread_continue(bounds, sizeof(double) * 6);
}
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           find_source_domain_action, find_source_domain_handler,
           HPX_POINTER, HPX_INT);


int find_target_domain_handler(Target *targets, int n_targets) {
  //NOTE: SMP assumptions

  //These are the three low bounds, followed by the three high bounds
  double bounds[6]{1.0e34, 1.0e34, 1.0e34, -1.0e34, -1.0e34, -1.0e34};

  for (int i = 0; i < n_targets; ++i) {
    Point p = targets[i].position();
    if (p.x() < bounds[0]) bounds[0] = p.x();
    if (p.x() > bounds[3]) bounds[3] = p.x();
    if (p.y() < bounds[1]) bounds[1] = p.y();
    if (p.y() > bounds[4]) bounds[4] = p.y();
    if (p.z() < bounds[2]) bounds[2] = p.z();
    if (p.z() > bounds[5]) bounds[5] = p.z();
  }

  hpx_thread_continue(bounds, sizeof(double) * 6);
}
HPX_ACTION(HPX_DEFAULT, HPX_PINNED,
           find_target_domain_action, find_target_domain_handler,
           HPX_POINTER, HPX_INT);


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

int evaluate_handler(EvaluateParams *parms, size_t total_size) {
  //copy user data into internal data
  hpx_addr_t source_packed = hpx_lco_future_new(sizeof(PackDataResult));
  assert(source_packed != HPX_NULL);
  hpx_addr_t target_packed = hpx_lco_future_new(sizeof(PackDataResult));
  assert(target_packed != HPX_NULL);

  hpx_call(parms->sources, pack_sources_action, source_packed,
           &parms->sources, &parms->spos_offset, &parms->q_offset);
  hpx_call(parms->targets, pack_targets_action, target_packed,
           &parms->targets, &parms->tpos_offset);

  hpx_addr_t source_bounds = hpx_lco_future_new(sizeof(double) * 6);
  assert(source_bounds != HPX_NULL);
  hpx_addr_t target_bounds = hpx_lco_future_new(sizeof(double) * 6);
  assert(target_bounds != HPX_NULL);

  hpx_call_when(source_packed, parms->sources, find_source_domain_action,
                source_bounds, &source_packed);
  hpx_call_when(target_packed, parms->targets, find_target_domain_action,
                target_bounds, &target_packed);

  //create our method and expansion from the parameters
  MethodSerial *method_serial = static_cast<MethodSerial *>(parms->data);
  hpx_addr_t method =
      hpx_gas_alloc_local_at_sync(1, parms->method_size, 0, HPX_HERE);
  hpx_gas_memput_rsync(method, method_serial, parms->method_size);

  ExpansionSerial *expansion_serial =
      static_cast<ExpansionSerial *>(parms->data + parms->method_size);
  hpx_addr_t expansion =
      hpx_gas_alloc_local_at_sync(1, parms->expansion_size, 0, HPX_HERE);
  hpx_gas_memput_rsync(expansion, expansion_serial, parms->expansion_size);

  //collect results of actions
  PackDataResult res{};
  hpx_lco_get(source_packed, sizeof(res), &res);
  hpx_lco_delete_sync(source_packed);
  int n_sources = res.count;
  hpx_addr_t source_parts = res.packed;

  hpx_lco_get(target_packed, sizeof(res), &res);
  hpx_lco_delete_sync(target_packed);
  int n_targets = res.count;
  hpx_addr_t target_parts = res.packed;

  DomainGeometry root_vol = cubify_domain(source_bounds, target_bounds);
  hpx_lco_delete_sync(source_bounds);
  hpx_lco_delete_sync(target_bounds);

  //build trees/do work
  SourceNode source_root{root_vol, 0, 0, 0, 0, method, nullptr};
  hpx_addr_t partitiondone =
      source_root.partition(source_parts, n_sources, parms->refinement_limit,
                            expansion);
  TargetNode target_root{root_vol, 0, 0, 0, 0, method, nullptr};
  hpx_lco_wait(partitiondone);
  hpx_lco_delete_sync(partitiondone);
  target_root.partition(target_parts, n_targets, parms->refinement_limit,
                        expansion, 0, std::vector<SourceNode>{});

  //copy results back into user data
  hpx_addr_t unpackdone = hpx_lco_future_new(0);
  assert(unpackdone != HPX_NULL);
  hpx_call(target_parts, unpack_targets_action, unpackdone,
           &target_parts, &parms->targets, &parms->phi_offset);

  //clean up
  source_root.destroy();
  target_root.destroy();
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
