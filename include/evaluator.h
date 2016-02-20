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


#ifndef __DASHMM_EVALUATOR_H__
#define __DASHMM_EVALUATOR_H__


#include "include/particle.h"


namespace dashmm {


template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename, typename> class Method>
class Evaluator {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Method<Source, Target, Expansion>;

  // Related type aliases to avoid clutter
  using sourceref_t = SourceRef<Source>;
  using targetref_t = TargetRef<Target>;
  using targetlco_t = TargetLCO<Source, Target, Expansion, Method>;
  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, Method>;
  using sourcenode_t = SourceNode<Source, Target, Expansion, Method>;
  using targetnode_t = TargetNode<Source, Target, Expansion, Method>;

  Evaluator() {
    // TODO: register all actions in the system

    // TargetLCO related
    HPX_REGISTER_ACTION(HPX_FUNCTION, 0,
                        targetlco_t::init_, targetlco_t::init_handler,
                        HPX_POINTER, HPX_SIZE_T, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_FUNCTION, 0,
                        targetlco_t::operation_,
                        targetlco_t::operation_handler,
                        HPX_POINTER, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_FUNCTION, 0,
                        targetlco_t::predicate_,
                        targetlco_t::predicate_handler,
                        HPX_POINTER, HPX_SIZE_T);

    // ExpansionLCO related
    HPX_REGISTER_ACTION(HPX_FUNCTION, 0,
                        targetlco_t::init_, targetlco_t::init_handler,
                        HPX_POINTER, HPX_SIZE_T, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_FUNCTION, 0,
                        expansionlco_t::operation_,
                        expansionlco_t::operation_handler,
                        HPX_POINTER, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_FUNCTION, 0,
                        expansionlco_t::predicate_,
                        expansionlco_t::predicate_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,
                        expansionlco_t::s_to_m_,
                        expansionlco_t::s_to_m_handler,
                        HPX_POINTER, HPX_INT, HPX_DOUBLE, HPX_DOUBLE,
                        HPX_DOUBLE, HPX_DOUBLE, HPX_ADDR, HPX_INT, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,
                        expansionlco_t::s_to_l_,
                        expansionlco_t::_s_to_l_handler,
                        HPX_POINTER, HPX_INT, HPX_DOUBLE, HPX_DOUBLE,
                        HPX_DOUBLE, HPX_DOUBLE, HPX_ADDR, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, 0,
                        expansionlco_t::m_to_m_,
                        expansionlco_t::m_to_m_handler,
                        HPX_ADDR, HPX_INT, HPX_INT, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        expansionlco_t::m_to_l_,
                        expansionlco_t::m_to_l_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, 0,
                        expansionlco_t::l_to_l_,
                        expansionlco_t::l_to_l_handler,
                        HPX_ADDR, HPX_INT, HPX_INT, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, 0,
                        expansionlco_t::m_to_t_,
                        expansionlco_t::m_to_t_handler,
                        HPX_INT, HPX_DOUBLE, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, 0,
                        expansionlco_t::l_to_t_,
                        expansionlco_t::l_to_t_handler,
                        HPX_INT, HPX_DOUBLE, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,
                        expansionlco_t::s_to_t_,
                        expansionlco_t::s_to_t_handler,
                        HPX_POINTER, HPX_INT, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, 0,
                        expansionlco_t::add_,
                        expansionlco_t::add_handler,
                        HPX_ADDR, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        expansionlco_t::create_from_expansion_,
                        expansionlco_t::create_from_expansion_handler,
                        HPX_POINTER, HPX_SIZE_T);

    // Source Node related
    HPX_REGISTER_ACTION(HPX_DEFAULT, 0,
                        sourcenode_t::self_delete_,
                        sourcenode_t::self_delete_handler,
                        HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, 0,
                        sourcenode_t::node_delete_,
                        sourcenode_t::node_delete_handler,
                        HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,
                        sourcenode_t::child_done_,
                        sourcenode_t::child_done_handler,
                        HPX_POINTER, HPX_ADDR, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED | HPX_MARSHALLED,
                        sourcenode_t::partition_,
                        sourcenode_t::partition_handler,
                        HPX_POINTER, HPX_POINTER, HPX_SIZE_T);

    // Target Node related
    HPX_REGISTER_ACTION(HPX_DEFAULT, 0,
                        targetnode_t::self_delete_,
                        targetnode_t::self_delete_handler,
                        HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, 0,
                        targetnode_t::node_delete_,
                        targetnode_t::node_delete_handler,
                        HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED | HPX_MARSHALLED,
                        targetnode_t::partition_,
                        targetnode_t::partition_handler,
                        HPX_POINTER, HPX_POINTER, HPX_SIZE_T);

    // Actions for the evaluation
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        evaluate_, evaluate_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, 0,
                        source_bounds_, source_bounds_handler,
                        HPX_ADDR, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, 0,
                        target_bounds_, target_bounds_handler,
                        HPX_ADDR, HPX_SIZE_T);
  }

  // TODO decide if we just make this operator()
  // NOTE: the arrays are marked const, even though the data will be sorted
  // and modified in both cases (just a sort of sources). The array object is
  // a reference. And so the reference passed in will not change. Only the
  // data that it is referring to will change.
  ReturnCode evaluate(const Array<source_t> sources,
                      const Array<target_t> targets,
                      int refinement_limit, const method_t &method,
                      const expansion_t &expansion) {
    // pack the arguments and call the action
    EvaluateParams args{ };
    args.sources = sources;
    args.targets = targets;
    args.refinement_limit = refinement_limit;
    args.method = method;
    args.n_digits = expansion.accuracy();

    if (HPX_SUCCESS != hpx_run(&evaluate_, &args, sizeof(args))) {
      return kRuntimeError;
    }

    return kSuccess;
  }

 private:
  hpx_action_t evaluate_;
  hpx_action_t source_bounds_;
  hpx_action_t target_bounds_;

  struct EvaluateParams {
    Array<source_t> sources;
    Array<target_t> targets;
    int refinement_limit;
    method_t method;
    int n_digits;
  };

  struct BoundsResult {
    Point low;
    Point high;
  };

  static int evaluate_handler(EvaluateParams *parms, size_t total_size) {
    // create source and target references
    // NOTE: These will need to be updated once we change the way array works
    // for distributed operation.
    ArrayMetaData srcmeta;
    hpx_gas_memget_sync(&srcmeta, parms->sources.data(), sizeof(srcmeta));
    sourceref_t sources{srcmeta.data, srcmeta.count, srcmeta.count};

    ArrayMetaData trgmeta;
    hpx_gas_memget_sync(&trgmeta, parms->targets.data(), sizeof(trgmeta));
    targetref_t targets{trgmeta.data, trgmeta.count, trgmeta.count};

    // find source and target bounds
    // TODO: Factor this bit
    hpx_addr_t srcbnd = hpx_lco_future_new(sizeof(BoundsResult));
    hpx_addr_t trgbnd = hpx_lco_future_new(sizeof(BoundsResult));
    assert(srcbnd != HPX_NULL && trgbnd != HPX_NULL);

    hpx_call(srcmeta.data, source_bounds_, srcbnd,
             &srcmeta.data, &srcmeta.count);
    hpx_call(trgmeta.data, target_bounds_, trgbnd,
             &trgmeta.data, &trgmeta.count);

    BoundsResult bounds{ };
    hpx_lco_get(srcbnd, sizeof(BoundsResult), &bounds);
    BoundsResult otherbounds{ };
    hpx_lco_get(trgbnd, sizeof(BoundsResult), &otherbounds);
    bounds.low.lower_bound(otherbounds.low);
    bounds.high.upper_bound(otherbounds.high);

    hpx_lco_delete_sync(srcbnd);
    hpx_lco_delete_sync(trgbnd);

    DomainGeometry domain = cubify_domain(bounds.low, bounds.high);

    // create source tree, wait for partitioning of source to finish,
    // partition target tree.
    sourcenode_t source_root{domain, Index{0, 0, 0}, parms->method, nullptr};
    hpx_addr_t partition_done =
      source_root.partition(sources, parms->refinement_limit, parms->n_digits);

    targetnode_t target_root{domain, Index{0, 0, 0}, parms->method, nullptr};
    hpx_lco_wait(partition_done);
    hpx_lco_delete_sync(partition_done);

    bool same_sandt = parms->sources.data() == parms->targets.data();
    target_root.partition(targets, parms->refinement_limint, parms->n_digits,
                          0, same_sandt,
                          std::vector<sourcenode_t>{source_root});

    // clean up
    source_root.destroy();
    target_root.destroy();

    // return
    hpx_exit(HPX_SUCCESS);
  }

  static int source_bounds_handler(hpx_addr_t data, size_t count) {
    source_t *user{nullptr};
    assert(hpx_gas_tru_pin(data, (void **)&user));

    BoundsResult retval{Point{1.0e34, 1.0e34, 1.0e34},
                        Point{-1.0e34, -1.0e34, -1.0e34}};

    for (size_t i = 0; i < count; ++i) {
      retval.low.lower_bound(user[i].position);
      retval.high.upper_bound(user[i].position);
    }

    return HPX_THREAD_CONTINUE(retval);
  }

  static int target_bounds_handler(hpx_addr_t data, size_t count) {
    target_t *user{nullptr};
    assert(hpx_gas_tru_pin(data, (void **)&user));

    BoundsResult retval{Point{1.0e34, 1.0e34, 1.0e34},
                        Point{-1.0e34, -1.0e34, -1.0e34}};

    for (size_t i = 0; i < count; ++i) {
      retval.low.lower_bound(user[i].position);
      retval.high.upper_bound(user[i].position);
    }

    return HPX_THREAD_CONTINUE(retval);
  }

  static int cubify_domain(Point low, Point high) {
    Point center = point_add(low.scale(0.5), high.scale(0.5));
    Point sizes = point_sub(high, low);
    double size = sizes.x() > sizes.y() ? sizes.x() : sizes.y();
    size = size > sizes.z() ? size : sizes.z();
    size *= 0.5001;
    Point offset{-size, -size, -size};

    return DomainGeometry{point_add(center, offset), size * 2.0};
  }
};


} //namespace dashmm


#endif // __DASHMM_EVALUATOR_H__
