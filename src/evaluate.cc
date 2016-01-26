// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  This software may be modified and distributed under the terms of the BSD
//  license. See the LICENSE file for details.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================


/// \file src/evaluate.cc
/// \brief Implementation of dashmm::evaluate()


#include <cstring>

#include <memory>

#include <hpx/hpx.h>

#include "include/array.h"
#include "include/expansion.h"
#include "include/method.h"
#include "include/methodref.h"
#include "include/node.h"
#include "include/types.h"


namespace dashmm {


/////////////////////////////////////////////////////////////////////
// Utility routines
/////////////////////////////////////////////////////////////////////


/// Cubifies the given domain boundaries
///
/// Given LCO addresses that hold the source and target domain boundaries,
/// this routine will compute a cubical volume that encompasses all source
/// and target points. The computed volume will extend slightly past the minmal
/// size that would contain the points.
///
/// \param s_bounds - an array that holds six double values that contain
///                        the domain bounds for the source points.
/// \param t_bounds - an array that holds six double values that contain
///                        the domain bounds for the target points.
///
/// \returns - a DomainGeometry that encompasses all sources and targets
DomainGeometry cubify_domain(double s_bounds[6],
                             double t_bounds[6]) {
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
  low = point_add(center, offset);
  size *= 2.0;

  return DomainGeometry{low, size};
}


/////////////////////////////////////////////////////////////////////
// Actions
/////////////////////////////////////////////////////////////////////


struct PackDataResult {
  hpx_addr_t packed;
  int count;
  double bounds[6]; 
};

/// Action that packs the needed source data and finds the bounding domain
///
/// To allow for the user's use of GAS to store more information with their
/// source records, DASHMM packs the required data into an internal structure.
/// This allows DASHMM to redistribute and sort the source data for performance.
///
/// This action continues the global address of the packed data and the count
/// of the packed records.
///
/// \param user_data - global address of the DASHMM array storing the user
///                    source data.
/// \param pos_offset - the offset in the source records to the three doubles
///                    giving the position of the source.
/// \param q_offset - the offset in the source records to the charge of the
///                   source
int pack_sources_handler(hpx_addr_t user_data, int pos_offset, int q_offset) {
  //NOTE: SMP assumptions all over this function
  ArrayMetaData *meta{nullptr};
  assert(hpx_gas_try_pin(user_data, (void **)&meta));

  char *user{nullptr};
  assert(hpx_gas_try_pin(meta->data, (void **)&user));

  PackDataResult retval{}; 
  retval.packed = hpx_gas_alloc_local_at_sync(1, meta->count * sizeof(Source), 
                                              0, HPX_HERE); 
  retval.count = meta->count; 

  double bounds[6]{1.0e34, 1.0e34, 1.0e34, -1.0e34, -1.0e34, -1.0e34}; 

  if (retval.packed != HPX_NULL) {
    Source *sources{nullptr}; 
    assert(hpx_gas_try_pin(retval.packed, (void **)&sources));

    for (size_t i = 0; i < meta->count; ++i) {
      void *pos_base = static_cast<void *>(&user[i * meta->size] + pos_offset);
      double *pos = static_cast<double *>(pos_base);
      void *q_base = static_cast<void *>(&user[i * meta->size] + q_offset);
      double *q = static_cast<double *>(q_base);
      sources[i].position = Point{pos[0], pos[1], pos[2]};
      sources[i].charge = *q;
      
      if (pos[0] < bounds[0]) bounds[0] = pos[0]; 
      if (pos[0] > bounds[3]) bounds[3] = pos[0]; 
      if (pos[1] < bounds[1]) bounds[1] = pos[1]; 
      if (pos[1] > bounds[4]) bounds[4] = pos[1]; 
      if (pos[2] < bounds[2]) bounds[2] = pos[2]; 
      if (pos[2] > bounds[5]) bounds[5] = pos[2]; 
    }
    hpx_gas_unpin(retval.packed);
  }

  for (size_t i = 0; i < 6; ++i) 
    retval.bounds[i] = bounds[i]; 

  hpx_gas_unpin(meta->data);
  hpx_gas_unpin(user_data);

  return HPX_THREAD_CONTINUE(retval);
}
HPX_ACTION(HPX_DEFAULT, 0,
           pack_sources_action, pack_sources_handler,
           HPX_ADDR, HPX_INT, HPX_INT);


/// Action that packs the needed target data and finds the bounding domain
///
/// To allow for the user's use of GAS to store more information with their
/// target records, DASHMM packs the required data into an internal structure.
/// This allows DASHMM to redistribute and sort the target data for performance.
///
/// This action continues the global address of the packed data and the count
/// of the packed records.
///
/// \param user_data - global address of the DASHMM array storing the user
///                    target data.
/// \param pos_offset - the offset in the target records to the three doubles
///                    giving the position of the target.
int pack_targets_handler(hpx_addr_t user_data, int pos_offset) {
  //NOTE: SMP assumptions
  ArrayMetaData *meta{nullptr};
  assert(hpx_gas_try_pin(user_data, (void **)&meta));

  char *user{nullptr};
  assert(hpx_gas_try_pin(meta->data, (void **)&user));

  PackDataResult retval{};
  retval.packed = hpx_gas_alloc_local_at_sync(1, meta->count * sizeof(Target), 
                                              0, HPX_HERE);
  retval.count = meta->count;

  double bounds[6]{1.0e34, 1.0e34, 1.0e34, -1.0e34, -1.0e34, -1.0e34}; 

  if (retval.packed != HPX_NULL) {
    Target *targets{nullptr};
    assert(hpx_gas_try_pin(retval.packed, (void **)&targets));

    for (size_t i = 0; i < meta->count; ++i) {
      void *pos_base = static_cast<void *>(&user[i * meta->size] + pos_offset);
      double *pos = static_cast<double *>(pos_base);
      targets[i].position = Point{pos[0], pos[1], pos[2]};
      targets[i].index = i;
      targets[i].phi = dcomplex_t{0.0};

      if (pos[0] < bounds[0]) bounds[0] = pos[0]; 
      if (pos[0] > bounds[3]) bounds[3] = pos[0]; 
      if (pos[1] < bounds[1]) bounds[1] = pos[1]; 
      if (pos[1] > bounds[4]) bounds[4] = pos[1]; 
      if (pos[2] < bounds[2]) bounds[2] = pos[2]; 
      if (pos[2] > bounds[5]) bounds[5] = pos[2]; 
    }
    hpx_gas_unpin(retval.packed);
  }

  for (size_t i = 0; i < 6; ++i) 
    retval.bounds[i] = bounds[i]; 

  hpx_gas_unpin(meta->data);
  hpx_gas_unpin(user_data);

  return HPX_THREAD_CONTINUE(retval);
}
HPX_ACTION(HPX_DEFAULT, 0,
           pack_targets_action, pack_targets_handler,
           HPX_ADDR, HPX_INT);


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

/// Action accomplishing the evaluation
///
/// This action takes a marshalled argument containing an EvaluateParams
/// structure.
///
/// This action continues nothing.
///
/// \param parms - the parameters
/// \param total_size - the size of the parameters in bytes
int evaluate_handler(EvaluateParams *parms, size_t total_size) {
  //copy user data into internal data
  hpx_addr_t source_packed = hpx_lco_future_new(sizeof(PackDataResult)); 
  hpx_addr_t target_packed = hpx_lco_future_new(sizeof(PackDataResult)); 
  assert(source_packed != HPX_NULL && target_packed != HPX_NULL); 

  hpx_call(parms->sources, pack_sources_action, source_packed, 
           &parms->sources, &parms->spos_offset, &parms->q_offset); 
  hpx_call(parms->targets, pack_targets_action, target_packed, 
           &parms->targets, &parms->tpos_offset); 

  //create our method and expansion from the parameters
  MethodSerial *method_serial = reinterpret_cast<MethodSerial *>(parms->data);
  auto local_method = create_method(method_serial->type, method_serial);
  MethodRef method = globalize_method(local_method, HPX_HERE);
  delete local_method;

  char *expansion_base = parms->data + parms->method_size;
  int *type = reinterpret_cast<int *>(expansion_base + sizeof(int));
  int *n_digits = reinterpret_cast<int *>(expansion_base + sizeof(int) * 2);
  char *local_copy = static_cast<char *>(malloc(parms->expansion_size));
  memcpy(local_copy, expansion_base, parms->expansion_size);
  auto local_expansion =
    interpret_expansion(*type, local_copy, parms->expansion_size, *n_digits);
  ExpansionRef expansion =
    globalize_expansion(std::move(local_expansion), HPX_HERE);
  expansion.finalize();

  // Collect result of actions
  PackDataResult s_res{}, t_res{}; 

  hpx_lco_get(source_packed, sizeof(s_res), &s_res); 
  hpx_lco_delete_sync(source_packed); 
  SourceRef sources{s_res.packed, s_res.count, s_res.count}; 

  hpx_lco_get(target_packed, sizeof(t_res), &t_res); 
  hpx_addr_t target_data = t_res.packed; 
  int target_count = t_res.count; 
  hpx_lco_delete_sync(target_packed); 

  DomainGeometry root_vol = cubify_domain(s_res.bounds, t_res.bounds);

  //build trees/do work - NOTE the awkwardness with source reference... This
  // really ought to be improved.
  SourceNode source_root{root_vol, Index{0, 0, 0, 0}, method.data(), nullptr};
  hpx_addr_t partitiondone =
    source_root.partition(sources, parms->refinement_limit,
                          expansion.type(), expansion.data(),
                          expansion.accuracy());
  
  TargetNode target_root{root_vol, Index{0, 0, 0, 0}, method.data(), nullptr};
  hpx_lco_wait(partitiondone);
  hpx_lco_delete_sync(partitiondone);
  bool same_sandt = (parms->sources == parms->targets);
  target_root.partition(target_data, target_count, parms->refinement_limit,
                        expansion, 0, same_sandt,
                        std::vector<SourceNode>{source_root});
  hpx_gas_free_sync(target_data);


  //copy results back into user data
  target_root.collect_results(parms->targets, parms->phi_offset);

  //clean up
  source_root.destroy();
  target_root.destroy();
  expansion.destroy();
  method.destroy();
  sources.destroy();

  //return
  hpx_exit(HPX_SUCCESS);
}
HPX_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
           evaluate_action, evaluate_handler,
           HPX_POINTER, HPX_SIZE_T);


/////////////////////////////////////////////////////////////////////
// Interface
/////////////////////////////////////////////////////////////////////


ReturnCode evaluate(ObjectHandle sources, int spos_offset, int q_offset,
                    ObjectHandle targets, int tpos_offset, int phi_offset,
                    int refinement_limit,
                    std::unique_ptr<Method> method,
                    std::unique_ptr<Expansion> expansion) {
  if (!method->compatible_with(expansion.get())) {
    return kIncompatible;
  }

  //Marshall our arguments
  size_t method_size = method->bytes();
  MethodSerial *method_serial = method->release();
  size_t expansion_size = expansion->bytes();
  char *expansion_serial = static_cast<char *>(expansion->release());
  size_t total_size = method_size + expansion_size + sizeof(EvaluateParams);
  EvaluateParams *args = static_cast<EvaluateParams *>(malloc(total_size));
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
  memcpy(args->data, method_serial, method_size);
  memcpy(args->data + method_size, expansion_serial, expansion_size);

  if (HPX_SUCCESS != hpx_run(&evaluate_action, args, total_size)) {
    return kRuntimeError;
  }

  free(method_serial);
  free(expansion_serial);
  free(args);

  return kSuccess;
}


} // namespace dashmm
