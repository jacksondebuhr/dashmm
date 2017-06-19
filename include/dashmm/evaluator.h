// =============================================================================
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2017, Trustees of Indiana University,
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


/// \file
/// \brief Definition of DASHMM Evaluator object


#include <iostream>
#include <memory>

#include <hpx/hpx.h>
#include <libhpx/libhpx.h>

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
/// Evaluator will cause a fatal error. In the future, this restriction may
/// be enforced by the library.
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
  using targetlco_t = TargetLCO<Source, Target, Expansion, Method>;
  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, Method>;
  using sourcenode_t = Node<Source>;
  using targetnode_t = Node<Target>;
  using dualtree_t = DualTree<Source, Target, Expansion, Method>;
  using distropolicy_t = typename method_t::distropolicy_t;

  /// The constuctor takes care of all action registration that DASHMM needs
  /// for one particular combination of Source, Target, Expansion and Method.
  /// Much of the registration occurs via Registrar objects, of which Evaluator
  /// had a number. Finally, a few Evaluator specific actions are registered
  /// in this constructor.
  Evaluator() : tlcoreg_{}, elcoreg_{}, snodereg_{}, tnodereg_{},
                streereg_{}, ttreereg_{}, dtreereg_{} {
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        create_tree_, create_tree_handler,
                        HPX_ADDR, HPX_ADDR, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        create_DAG_, create_DAG_handler,
                        HPX_ADDR, HPX_INT, HPX_POINTER, HPX_INT,
                        HPX_POINTER, HPX_POINTER);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        execute_DAG_, execute_DAG_handler,
                        HPX_ADDR, HPX_POINTER);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        reset_DAG_, reset_DAG_handler,
                        HPX_POINTER);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        destroy_DAG_, destroy_DAG_handler,
                        HPX_ADDR. HPX_POINTER);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        destroy_tree_, destroy_tree_handler,
                        HPX_ADDR);
  }

  TreeHandle create_tree(const Array<source_t> &sources,
                         const Array<source_t> &targets,
                         int refinement_limit) {
    hpx_addr_t sources_addr{sources.data()};
    hpx_addr_t targets_addr{targets.data()};
    hpx_addr_t rwaddr{HPX_NULL};
    hpx_run(&create_tree_, &rwaddr, &sources_addr, &targets_addr,
            &refinement_limit);

    return rwaddr;
  }

  std::unique_ptr<DAG> create_DAG(TreeHandle tree,
                                  const method_t *method,
                                  int n_digits,
                                  const std::vector<double> *kernel_params,
                                  distropolicy_t distro = distropolicy_t{}) {
    int n_params{kernel_params.size()};
    const distropolicy_t *distro_ptr{&distro};
    DAG *dag{nullptr};
    hpx_run_spmd(&create_DAG_, &dag, &tree, &n_digits, &kernel_params,
                 &n_params, &distro_ptr, &method);
    return std::unique_ptr<DAG>{dag};
  }

  ReturnCode execute_DAG(TreeHandle tree, DAG *dag) {
    if (HPX_SUCCESS == hpx_run_spmd(&execute_DAG_, nullptr, &tree, &dag)) {
      return kSuccess;
    } else {
      return kRuntimeError;
    }
  }

  ReturnCode reset_DAG(DAG *dag) {
    if (HPX_SUCCESS == hpx_run_spmd(&reset_DAG_, nullptr, &dag)) {
      return kSuccess;
    } else {
      return kRuntimeError;
    }
  }

  // TODO: The need to include the tree here is a sign that something is off.
  // In particular, the destroy LCOs thing should not need to know about the
  // tree. I just have that as a member of tree because that file got out
  // of hand.
  ReturnCode destroy_DAG(TreeHandle tree, std::unique_ptr<DAG> dag) {
    DAG *ptr = dag.release();
    if (HPX_SUCCESS == hpx_run_spmd(&destroy_DAG_, nullptr, &tree, &ptr)) {
      return kSuccess;
    } else {
      return kRuntimeError;
    }
  }

  ReturnCode destroy_tree(TreeHandle rwaddr) {
    if (HPX_SUCCESS == hpx_run(&destroy_tree_, nullptr, &rwaddr)) {
      return kSuccess;
    } else {
      return kRuntimeError;
    }
  }

  /// Perform a multipole moment evaluation
  ///
  /// Given source, targets, a refinement_limit, a method, the accuracy
  /// parameter, the kernel parameters, and an optional distribution policy,
  /// this routine performs a multipole method evaluation. The source and
  /// target arrays will likely be sorted by DASHMM during evaluation, so it
  /// is important that users include some sort of identifying data in the
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
  /// The provided method is used as prototypes for any other copies of those
  /// objects that are required. Some methods have data associated with them.
  /// In general, expansions might need an accuracy parameter and some other
  /// parameters that control the potential. These are passed in as @p n_digits
  /// and @p kernelparams.
  ///
  /// Finally, @distro is an optional parameter if the user would like the
  /// provide a policy that is not default constructed.
  ///
  /// \param sources - a DASHMM Array of the source points
  /// \param targets - a DASHMM Array of the target points
  /// \param refinement_limit - the domain refinement limit
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
                      int refinement_limit, const method_t *method,
                      int n_digits, const std::vector<double> *kernelparams,
                      distropolicy_t distro = distropolicy_t{}) {
    TreeHandle tree = create_tree(sources, targets, refinement_limit);
    auto dag = create_DAG(tree, method, n_digits, kernelparams, distro);
    if (kRuntimeError == execute_DAG(tree, dag.get())) {
      // TODO: what do we do about that? How to clean up?
      // These are a real problem.
      return kRuntimeError;
    }
    if (kRuntimeError == destroy_DAG(tree, std::move(dag))) {
      // TODO: what do we do about that?
      return kRuntimeError;
    }
    if (kRuntimeError == destroy_tree(tree)) {
      // TODO: what do we do about that?
      return kRuntimeError;
    }

    return kSuccess;
  }

 private:
  /// Registrars that will be needed for evaluate to operate.
  TargetLCORegistrar<Source, Target, Expansion, Method> tlcoreg_;
  ExpansionLCORegistrar<Source, Target, Expansion, Method> elcoreg_;
  NodeRegistrar<Source> snodereg_;
  NodeRegistrar<Target> tnodereg_;
  TreeRegistrar<Source, Target, Source, Expansion, Method> streereg_;
  TreeRegistrar<Source, Target, Target, Expansion, Method> ttreereg_;
  DualTreeRegistrar<Source, Target, Expansion, Method> dtreereg_;

  // The actions for evaluate
  static hpx_action_t create_tree_;
  static hpx_action_t create_DAG_;
  static hpx_action_t execute_DAG_;
  static hpx_action_t reset_DAG_;
  static hpx_action_t destroy_DAG_;
  static hpx_action_t destroy_tree_;

  static int create_tree_handler(hpx_addr_t sources_addr,
                                 hpx_addr_t targets_addr,
                                 int refinement_limit) {
    Array<source_t> sources{sources_addr};
    Array<target_t> targets{targets_addr};

    hpx_time_t creation_begin = hpx_time_now();
    RankWise<dualtree_t> global_tree =
    dualtree_t::create(refinement_limit, sources, targets);

    hpx_addr_t partitiondone =
        dualtree_t::partition(global_tree, sources, targets);
    hpx_lco_wait(partitiondone);
    hpx_lco_delete_sync(partitiondone);
    hpx_time_t creation_end = hpx_time_now();
    double creation_deltat = hpx_time_diff_us(creation_begin, creation_end);
    fprintf(stdout, "Evaluation: tree creation %lg [us]\n", creation_deltat);

    hpx_addr_t rwaddr = global_tree.data();
    hpx_exit(sizeof(hpx_addr_t), &rwaddr);
  }

  static int create_DAG_handler(hpx_addr_t rwaddr, int n_digits,
                                const std::vector<double> *kernel_params,
                                int n_params,
                                const distropolicy_t *distro_ptr,
                                const method_t *method_ptr) {
    distropolicy_t distro{*distro_ptr};
    method_t method{*method_ptr};

    // Generate the table
    RankWise<dualtree_t> global_tree{rwaddr};
    auto tree = global_tree.here();
    double domain_size = tree->domain()->size();
    expansion_t::update_table(n_digits, domain_size, *kernel_params);
    tree->set_method(method);

    // This creates and distributes the explicit DAG
    hpx_time_t distribute_begin = hpx_time_now();
    DAG *dag = tree->create_DAG();
    distro.compute_distribution(*dag);
    hpx_time_t distribute_end = hpx_time_now();
    double distribute_deltat = hpx_time_diff_us(distribute_begin,
                                                distribute_end);
    fprintf(stdout, "Evaluate: DAG creation and distribution: %lg [us]\n",
            distribute_deltat);

    // This allocates the implicit DAG
    hpx_time_t allocate_begin = hpx_time_now();
    tree->create_expansions_from_DAG(rwaddr);
    hpx_time_t allocate_end = hpx_time_now();
    double allocate_deltat = hpx_time_diff_us(allocate_begin, allocate_end);
    fprintf(stdout, "Evaluate: LCO allocation: %lg [us]\n", allocate_deltat);

    //return a pointer to the DAG at this rank
    hpx_exit(sizeof(DAG *), &dag);
  }

  static int execute_DAG_handler(hpx_addr_t rwaddr, DAG *dag,
                                 hpx_addr_t middone) {
    RankWise<dualtree_t> global_tree{rwaddr};
    auto tree = global_tree.here();

#ifdef DASHMM_INSTRUMENTATION
    libhpx_inst_phase_begin();
    // This is repeated, because the previous call is an atomic operation with
    // a relaxed memory model. That means the first of these only sometimes
    // will emit an event. In (simple and likely incomplete) testing, two
    // is sufficient to always emit at least one.
    EVENT_TRACE_DASHMM_ZEROREF();
    EVENT_TRACE_DASHMM_ZEROREF();
#endif
    hpx_time_t evaluate_begin = hpx_time_now();
    tree->start_DAG_evaluation(global_tree, dag);

    hpx_addr_t heredone = tree->setup_termination_detection(dag);
    hpx_lco_wait(heredone);

    hpx_time_t evaluate_end = hpx_time_now();
    double evaluate_deltat = hpx_time_diff_us(evaluate_begin, evaluate_end);
    fprintf(stdout, "Evaluate: DAG evaluation: %lg [us]\n", evaluate_deltat);

#ifdef DASHMM_INSTRUMENTATION
    libhpx_inst_phase_end();
#endif

    hpx_lco_delete_sync(heredone);
    hpx_lco_delete_sync(middone);
    hpx_exit(0, nullptr);
  }

  static int reset_DAG_handler(DAG *dag) {
    // TODO implement
    // basically a routine that will go ahead and reset any local LCO in
    // the ExpansionLCO or TargetLCO kinds

    hpx_exit(0, nullptr);
  }

  static int destroy_DAG_handler(hpx_addr_t rwaddr, DAG *dag) {
    RankWise<dualtree_t> global_tree{rwaddr};
    auto tree = global_tree.here();
    tree->destroy_DAG_LCOs(*dag);
    delete dag;
    hpx_exit(0, nullptr);
  }

  static int destroy_tree_handler(hpx_addr_t rwaddr) {
    RankWise<dualtree_t> global_tree{rwaddr};
    dualtree_t::destroy(global_tree);
    hpx_exit(0, nullptr);
  }
};


template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Evaluator<S, T, E, M>::create_tree_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Evaluator<S, T, E, M>::create_DAG_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Evaluator<S, T, E, M>::execute_DAG_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Evaluator<S, T, E, M>::reset_DAG_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Evaluator<S, T, E, M>::destroy_DAG_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Evaluator<S, T, E, M>::destroy_tree_ = HPX_ACTION_NULL;

} // namespace dashmm


#endif // __DASHMM_EVALUATOR_H__
