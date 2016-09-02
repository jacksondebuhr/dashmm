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
#include "dashmm/arrayref.h"
#include "dashmm/defaultpolicy.h"
#include "dashmm/domaingeometry.h"
#include "dashmm/expansionlco.h"
#include "dashmm/point.h"
#include "dashmm/rankwise.h"
#include "dashmm/registrar.h"
#include "dashmm/targetlco.h"
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
                    template <typename, typename> class,
                    typename> class Method,
          typename DistroPolicy = DefaultDistributionPolicy>
class Evaluator {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Method<Source, Target, Expansion, DistroPolicy>;
  using sourceref_t = ArrayRef<Source>;
  using targetref_t = ArrayRef<Target>;
  using targetlco_t = TargetLCO<Source, Target, Expansion, Method,
                                DistroPolicy>;
  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, Method,
                                      DistroPolicy>;
  using sourcenode_t = Node<Source, Target, Source, Expansion, Method,
                                DistroPolicy>;
  using targetnode_t = Node<Source, Target, Target, Expansion, Method,
                                DistroPolicy>;
  using dualtree_t = DualTree<Source, Target, Expansion, Method, DistroPolicy>;
  using distropolicy_t = DistroPolicy;

  /// The constuctor takes care of all action registration that DASHMM needs
  /// for one particular combination of Source, Target, Expansion and Method.
  Evaluator() : tlcoreg_{}, elocreg_{}, snodereg_{}, tnodereg_{},
                streereg_{}, ttreereg_{}, dtreereg_{} {
    // Actions for the evaluation
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        evaluate_, evaluate_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        evaluate_rank_local_, evaluate_rank_local_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        evaluate_cleanup_, evaluate_cleanup_handler,
                        HPX_ADDR, HPX_ADDR);
  }

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
  /// \param n_digits - the number of digits of accuracy required
  /// \param kernelparams - the parameters needed by the kernel
  /// \param distro - an instance of the distribution policy to use for this
  ///                 execution.
  ///
  /// \returns - kSuccess on success; kRuntimeError if there is an error with
  ///            the runtime.
  ReturnCode evaluate(const Array<source_t> &sources,
                      const Array<target_t> &targets,
                      int refinement_limit, const method_t &method,
                      int n_digits, const std::vector<double> &kernelparams,
                      const distropolicy_t &distro = distropolicy_t { }) {
    // pack the arguments and call the action
    size_t n_params = kernelparams.size();
    size_t total_size = sizeof(EvaluateParams) + n_params * sizeof(double);
    EvaluateParams *args =
        reinterpret_cast<EvaluateParams *>(new char[total_size]);
    args->sources = sources;
    args->targets = targets;
    args->refinement_limit = refinement_limit;
    args->method = method;
    args->n_digits = n_digits;
    args->distro = distro;
    for (size_t i = 0; i < n_params; ++i) {
      args->kernelparams[i] = kernelparams[i];
    }

    if (HPX_SUCCESS != hpx_run(&evaluate_, nullptr, args, total_size)) {
      return kRuntimeError;
    }

    return kSuccess;
  }

 private:
  /// Registrars that will be needed for evaluate to operate.
  TargetLCORegistrar<Source, Target, Expansion, Method,
                     DistroPolicy> tlcoreg_;
  ExpansionLCORegistrar<Source, Target, Expansion, Method,
                        DistroPolicy> elcoreg_;
  NodeRegistrar<Source> snodereg_;
  NodeRegistrar<Target> tnodereg_;
  TreeRegistrar<Source, Target, Source, Expansion, Method,
                DistroPolicy> streereg_;
  TreeRegistrar<Source, Target, Target, Expansion, Method,
                DistroPolicy> ttreereg_;
  DualTreeRegistrar<Source, Target, Expansion, Method, DistroPolicy> dtreereg_;

  // The actions for evaluate
  static hpx_action_t evaluate_;
  static hpx_action_t evaluate_rank_local_;
  static hpx_action_t evaluate_cleanup_;

  /// Parameters to evaluations
  struct EvaluateParams {
    Array<source_t> sources;
    Array<target_t> targets;
    size_t refinement_limit;
    method_t method;
    int n_digits;
    distropolicy_t distro;
    double kernelparams[];
  };

  /// The evaluation action implementation
  static int evaluate_handler(EvaluateParams *parms, size_t total_size) {
    RankWise<dualtree_t> global_tree =
        dualtree_t::create(parms->refinement_limit, parms->sources,
                           parms->targets);
    dualtree_t::partition(global_tree, sources, targets);

    // Allocate space for message buffer
    size_t bcast_size = total_size + sizeof(hpx_addr_t) * 2;
    char *data = new char[bcast_size];
    assert(data != nullptr);
    hpx_addr_t *preargs = reinterpret_cast<hpx_addr_t *>(data);
    EvaluateParams *postargs
        = reinterpret_cast<EvaluateParams *>(data + sizeof(hpx_addr_t) * 2);

    // Setup message
    preargs[0] = hpx_lco_and_new(hpx_get_num_ranks());
    assert(preargs[0] != HPX_NULL);
    preargs[1] = global_tree.data();
    memcpy(postargs, parms, total_size);

    hpx_bcast_lsync(evaluate_rank_local_, HPX_NULL, data, bcast_size);

    // set up dependent call on the broadcast to do evaluate cleanup
    hpx_call_when(alldone, HPX_HERE, evaluate_cleanup_,
                  &preargs[1], &preargs[0]);

    delete [] data;

    return HPX_SUCCESS;
  }

  // This is a marshalled action that is the broadcast target
  static int evaluate_rank_local_handler(char *data, size_t msg_size) {
    hpx_addr_t *preargs = reinterpret_cast<hpx_addr_t *>(data);
    EvaluateParams *parms
        = reinterpret_cast<EvaluateParams *>(data + sizeof(hpx_addr_t) * 2);

    // Generate the table - TODO: In principle, this can begin as soon as
    // the domain geometry is ready. However, starting it then might be a
    // difficult thing to arrange.
    RankWise<dualtree_t> global_tree{preargs[1]};
    auto tree = global_tree.here();
    double domain_size = tree->domain_.size();
    size_t n_params = (total_size - sizeof(EvaluateParams)) / sizeof(double);
    std::vector<double> kernel_params(parms->kernelparams,
                                      &parms->kernelparams[n_params]);
    expansion_t::update_table(parms->n_digits, domain_size, kernel_params);

    // Get ready to evaluate
    DAG dag = tree->create_DAG();
    parms->distro.compute_distribution(dag);
    tree->create_expansions_from_DAG();

    // NOTE: the previous has to finish for the following. So the previous
    // is a synchronous operation. The next three, however, are not. They all
    // get their work going when they come to it and then they return.

    tree->setup_edge_lists(dag);
    tree->start_DAG_evaluation();
    hpx_addr_t heredone = tree->setup_termination_detection(dag);

    // When the local work finished, delete the gate, and then set the
    // global completion LCO.
    hpx_call_when_with_continuation(heredone, heredone, hpx_lco_delete_action,
                                    &preargs[0], hpx_lco_set_action,
                                    nullptr, 0);

    return HPX_SUCCESS;
  }

  static int evaluate_cleanup_handler(hpx_addr_t rwaddr, hpx_addr_t alldone) {
    hpx_lco_delete_sync(alldone);

    // clean up tree and table
    RankWise<dualtree_t> global_tree{rwaddr};
    dualtree_t::destroy(global_tree);

    // Exit from this HPX epoch
    hpx_exit(0, nullptr);
  }
};


template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Evaluator<S, T, E, M, D>::evaluate_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Evaluator<S, T, E, M, D>::evaluate_rank_local_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Evaluator<S, T, E, M, D>::evaluate_cleanup_ = HPX_ACTION_NULL;


} // namespace dashmm


#endif // __DASHMM_EVALUATOR_H__
