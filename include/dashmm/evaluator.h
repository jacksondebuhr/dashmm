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


/// \file include/dashmm/evaluator.h
/// \brief Definition of DASHMM Evaluator object


#include <hpx/hpx.h>

#include "dashmm/array.h"
#include "dashmm/domaingeometry.h"
#include "dashmm/expansionlco.h"
#include "dashmm/point.h"
#include "dashmm/sourceref.h"
#include "dashmm/targetlco.h"
#include "dashmm/targetref.h"
#include "dashmm/tree.h"



namespace dashmm {


/// Evaluator object
///
/// This object is a central object in DASHMM evaluations. This object bears
/// two responsibilities: performing multipole method evaluations, and
/// registration of actions with the runtime system.
///
/// The interface for the object contains only one member, evaluate(), which
/// performs a multipole method evaluation.
///
/// In addition, the object's constructor performs its second duty. DASHMM is
/// a templated library. HPX-5 requires the address of functions that are to
/// be called as actions. Thus, we must somehow manage to register the specific
/// instances of the actions DASHMM requires once the template arguments are
/// specified. By creating an Evaluator object, which requires filling out the
/// template arguments, this object is able to instantiate all the templates
/// needed by the particular evaluation represented by this object.
///
/// A further note on the use of Evaluator. To register a particular set of
/// types with the runtime, an instance of this object must be created before
/// dashmm::init(). Evaluators created after dashmm::init() will cause a
/// fatal error.
///
/// Finally, only one instance of an evaluator for a given set of types
/// should be created. Multiple instances of a particular instance of an
/// Evaluator will cause a fatal error.
///
/// In the future, this restriction may be enforced by the library.
///
/// This is a template class over the Source, Targer, Expansion and Method
/// types. For a description of the requirements of these types, please see
/// the documentation.
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class> class Method>
class Evaluator {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Method<Source, Target, Expansion>;
  using sourceref_t = SourceRef<Source>;
  using targetref_t = TargetRef<Target>;
  using targetlco_t = TargetLCO<Source, Target, Expansion, Method>;
  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, Method>;
  using sourcenode_t = TreeNode<Source, Target, Source, Expansion, Method>;
  using targetnode_t = TreeNode<Source, Target, Target, Expansion, Method>;
  using tree_t = Tree<Source, Target, Expansion, Method>;

  /// The constuctor takes care of all action registration that DASHMM needs
  /// for one particular combination of Source, Target, Expansion and Method.
  ///
  /// For this to work, the related classes must mark Evaluator as a friend.
  Evaluator() {
    // TargetLCO related
    HPX_REGISTER_ACTION(HPX_FUNCTION, HPX_ATTR_NONE,
                        targetlco_t::init_, targetlco_t::init_handler,
                        HPX_POINTER, HPX_SIZE_T, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_FUNCTION, HPX_ATTR_NONE,
                        targetlco_t::operation_,
                        targetlco_t::operation_handler,
                        HPX_POINTER, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_FUNCTION, HPX_ATTR_NONE,
                        targetlco_t::predicate_,
                        targetlco_t::predicate_handler,
                        HPX_POINTER, HPX_SIZE_T);

    // ExpansionLCO related
    HPX_REGISTER_ACTION(HPX_FUNCTION, HPX_ATTR_NONE,
                        expansionlco_t::init_, expansionlco_t::init_handler,
                        HPX_POINTER, HPX_SIZE_T, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_FUNCTION, HPX_ATTR_NONE,
                        expansionlco_t::operation_,
                        expansionlco_t::operation_handler,
                        HPX_POINTER, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_FUNCTION, HPX_ATTR_NONE,
                        expansionlco_t::predicate_,
                        expansionlco_t::predicate_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,
                        expansionlco_t::s_to_m_,
                        expansionlco_t::s_to_m_handler,
                        HPX_POINTER, HPX_INT, HPX_DOUBLE, HPX_DOUBLE,
                        HPX_DOUBLE, HPX_DOUBLE, HPX_ADDR, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,
                        expansionlco_t::s_to_l_,
                        expansionlco_t::s_to_l_handler,
                        HPX_POINTER, HPX_INT, HPX_DOUBLE, HPX_DOUBLE,
                        HPX_DOUBLE, HPX_DOUBLE, HPX_ADDR, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        expansionlco_t::m_to_m_,
                        expansionlco_t::m_to_m_handler,
                        HPX_ADDR, HPX_INT, HPX_INT, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        expansionlco_t::m_to_l_,
                        expansionlco_t::m_to_l_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        expansionlco_t::l_to_l_,
                        expansionlco_t::l_to_l_handler,
                        HPX_ADDR, HPX_INT, HPX_INT, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        expansionlco_t::m_to_t_,
                        expansionlco_t::m_to_t_handler,
                        HPX_INT, HPX_DOUBLE, HPX_ADDR, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        expansionlco_t::l_to_t_,
                        expansionlco_t::l_to_t_handler,
                        HPX_INT, HPX_DOUBLE, HPX_ADDR, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,
                        expansionlco_t::s_to_t_,
                        expansionlco_t::s_to_t_handler,
                        HPX_POINTER, HPX_INT, HPX_ADDR, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        expansionlco_t::add_,
                        expansionlco_t::add_handler,
                        HPX_ADDR, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        expansionlco_t::create_from_expansion_,
                        expansionlco_t::create_from_expansion_handler,
                        HPX_POINTER, HPX_SIZE_T);

    // Tree related
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        tree_t::source_child_done_,
                        tree_t::source_child_done_handler,
                        HPX_POINTER, HPX_POINTER, HPX_ADDR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        tree_t::source_partition_,
                        tree_t::source_partition_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        tree_t::target_partition_,
                        tree_t::target_partition_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        tree_t::source_bounds_, tree_t::source_bounds_handler,
                        HPX_ADDR, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        tree_t::target_bounds_, tree_t::target_bounds_handler,
                        HPX_ADDR, HPX_SIZE_T);

    // Actions for the evaluation
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        evaluate_, evaluate_handler,
                        HPX_POINTER, HPX_SIZE_T);

  }


  // NOTE: the arrays are marked const, even though the data will be sorted
  // and modified in both cases (just a sort of sources). The array object is
  // a reference. And so the reference passed in will not change. Only the
  // data that it is referring to will change.

  /// Perform a multipole moment evaluation
  ///
  /// Given source, targets, a refinement_limit, a method and an expansion,
  /// this routine performs a multipole method evaluation. The source and
  /// target arrays will likely be sorted by DASHMM during evaluation, so it
  /// is imperative that users include some sort of identifying data in the
  /// given arrays.
  ///
  /// In addition to sorting the data, the records in the targets Array will
  /// be given the results of the evaluation. Other than these two changes,
  /// the arrays will be unchanged. NOTE: the sources and targets are passed
  /// as const even though the data they represent will be modified. Array
  /// objects are references to data in the global address space. It is the
  /// global data that will be modified, not where that data can be found. For
  /// this reason, we can pass these arguments as const.
  ///
  /// The refinement limit controls the partitioning of the domain. A given
  /// portion of the domain will be subdivided if there are more than the
  /// given number of sources or targets in that region.
  ///
  /// The provided method and expansion are used as prototypes for any other
  /// copies of those objects that are required. Some methods have data
  /// associated with them. Some expansions require an accuracy parameter to
  /// fully specify the expansion, and so this argument will supply that
  /// accuracy.
  ///
  /// \param sources - a DASHMM Array of the source points
  /// \param targets - a DASHMM Array of the target points
  /// \param refinement_limint - the domain refinement limit
  /// \param method - a prototype of the method to use.
  /// \param expansion - a prototpe of the expansion to use.
  ///
  /// \returns - kSuccess on success; kRuntimeError if there is an error with
  ///            the runtime.
  ReturnCode evaluate(const Array<source_t> &sources,
                      const Array<target_t> &targets,
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
  /// Action for evalutation
  static hpx_action_t evaluate_;
  /// Action for finding the source bounds
  static hpx_action_t source_bounds_;
  /// Action for finding the target bounds
  static hpx_action_t target_bounds_;

  /// Parameters to evaluations
  struct EvaluateParams {
    Array<source_t> sources;
    Array<target_t> targets;
    int refinement_limit;
    method_t method;
    int n_digits;
  };

  /// The result of finding bounds
  struct BoundsResult {
    Point low;
    Point high;
  };

  /// The evaluation action implementation
  static int evaluate_handler(EvaluateParams *parms, size_t total_size) {
    bool same_sandt = parms->sources.data() == parms->targets.data();

    // create source and target references
    // NOTE: These will need to be updated once we change the way array works
    // for distributed operation.
    ArrayMetaData srcmeta;
    hpx_gas_memget_sync(&srcmeta, parms->sources.data(), sizeof(srcmeta));
    sourceref_t sources{srcmeta.data, srcmeta.count, srcmeta.count};

    ArrayMetaData trgmeta;
    hpx_gas_memget_sync(&trgmeta, parms->targets.data(), sizeof(trgmeta));
    targetref_t targets{trgmeta.data, trgmeta.count, trgmeta.count};

    // Get source bounds
    hpx_addr_t srcbnd = hpx_lco_future_new(sizeof(BoundsResult));
    assert(srcbnd != HPX_NULL);
    hpx_call(srcmeta.data, source_bounds_, srcbnd,
             &srcmeta.data, &srcmeta.count);

    BoundsResult bounds{Point{0.0, 0.0, 0.0}, Point{0.0, 0.0, 0.0}};
    hpx_lco_get(srcbnd, sizeof(BoundsResult), &bounds);
    hpx_lco_delete_sync(srcbnd);

    // get target bounds - if targets != sources
    if (!same_sandt) {
      hpx_addr_t trgbnd = hpx_lco_future_new(sizeof(BoundsResult));
      assert(trgbnd != HPX_NULL);

      hpx_call(trgmeta.data, target_bounds_, trgbnd,
               &trgmeta.data, &trgmeta.count);

      BoundsResult otherbounds{Point{0.0, 0.0, 0.0}, Point{0.0, 0.0, 0.0}};
      hpx_lco_get(trgbnd, sizeof(BoundsResult), &otherbounds);
      bounds.low.lower_bound(otherbounds.low);
      bounds.high.upper_bound(otherbounds.high);

      hpx_lco_delete_sync(trgbnd);
    }

    // Get the overall domain
    DomainGeometry domain{bounds.low, bounds.high, 1.0002};

    // create source tree, wait for partitioning of source to finish,
    // partition target tree.
    sourcenode_t source_root{domain, Index{0, 0, 0, 0}, parms->method, nullptr};
    hpx_addr_t partition_done =
      source_root.partition(sources, parms->refinement_limit, parms->n_digits);

    targetnode_t target_root{domain, Index{0, 0, 0, 0}, parms->method, nullptr};
    hpx_lco_wait(partition_done);
    hpx_lco_delete_sync(partition_done);

    hpx_addr_t targetpartdone =
      target_root.partition(targets, parms->refinement_limit,
                            parms->n_digits, 0, same_sandt,
                            std::vector<sourcenode_t>{source_root});


    // deal with one pathological case:
    expansionlco_t srootexpand = source_root.expansion();
    hpx_lco_wait(srootexpand.data());

    // deal with another pathological case:
    hpx_lco_wait(targetpartdone);
    hpx_lco_delete_sync(targetpartdone);

    // clean up
    source_root.destroy();
    target_root.destroy();

    // return
    hpx_exit(HPX_SUCCESS);
  }
};


template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Evaluator<S, T, E, M>::evaluate_ = HPX_ACTION_NULL;


} // namespace dashmm


#endif // __DASHMM_EVALUATOR_H__
