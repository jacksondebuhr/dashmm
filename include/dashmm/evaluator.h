// =============================================================================
//  This file is part of:
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  DASHMM is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  DASHMM is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with DASHMM. If not, see <http://www.gnu.org/licenses/>.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================


#ifndef __DASHMM_EVALUATOR_H__
#define __DASHMM_EVALUATOR_H__


/// \file
/// \brief Definition of DASHMM Evaluator object


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
    // Actions for the evaluation
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        evaluate_, evaluate_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,
                        evaluate_rank_local_, evaluate_rank_local_handler,
                        HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,
                        evaluate_cleanup_, evaluate_cleanup_handler,
                        HPX_ADDR, HPX_ADDR, HPX_ADDR);
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
    args->rwaddr = HPX_NULL;
    args->alldone = HPX_NULL;
    args->middone = HPX_NULL;
    for (size_t i = 0; i < n_params; ++i) {
      args->kernelparams[i] = kernelparams[i];
    }

    if (HPX_SUCCESS != hpx_run(&evaluate_, nullptr, args, total_size)) {
      return kRuntimeError;
    }

    delete [] args;

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
    hpx_addr_t rwaddr;
    hpx_addr_t alldone;
    hpx_addr_t middone;
    double kernelparams[];
  };

  /// The evaluation action implementation
  ///
  /// This action is a diffusive action that starts the evaluation process.
  /// Here, all of the collective work is begun. In particular, the tree is
  /// constructed and a few synchronization LCOs are created. After this,
  /// a few dependent actions are set up for when the various operations are
  /// complete.
  ///
  /// NOTE: This action does not call hpx_exit().
  ///
  /// \param parms - the input arguments
  /// \param total_size - the size of the arguments
  ///
  /// \returns - HPX_SUCCESS
  static int evaluate_handler(EvaluateParams *parms, size_t total_size) {
    // BEGIN TREE CREATION
    //hpx_time_t creation_begin = hpx_time_now();
    RankWise<dualtree_t> global_tree =
        dualtree_t::create(parms->refinement_limit, parms->sources,
                           parms->targets);
    hpx_addr_t partitiondone =
        dualtree_t::partition(global_tree, parms->sources, parms->targets);
    hpx_lco_wait(partitiondone);
    hpx_lco_delete_sync(partitiondone);
    //hpx_time_t creation_end = hpx_time_now();
    //double creation_deltat = hpx_time_diff_us(creation_begin, creation_end);
    //fprintf(stdout, "Evaluation: tree creation %lg [us]\n", creation_deltat);
    // END TREE CREATION

    //auto local_tree = global_tree.here();
    //fprintf(stdout, "UnifGridDistrib:");
    //for (int i = 0; i < hpx_get_num_ranks(); ++i) {
    //  fprintf(stdout, " %d", local_tree->last(i) - local_tree->first(i) + 1);
    //}
    //fprintf(stdout, "\n");

    // Save tree global address in message
    parms->rwaddr = global_tree.data();
    // Save alldone LCO into message
    parms->alldone = hpx_lco_and_new(hpx_get_num_ranks());
    assert(parms->alldone != HPX_NULL);
    // The third is the expansion creation done LCO
    parms->middone = hpx_lco_and_new(hpx_get_num_ranks());
    assert(parms->middone != HPX_NULL);

    // Start the work everywhere
    hpx_bcast_lsync(evaluate_rank_local_, HPX_NULL, parms, total_size);

    // set up dependent call on the broadcast to do evaluate cleanup
    hpx_call_when(parms->alldone, HPX_HERE, evaluate_cleanup_, HPX_NULL,
                  &parms->rwaddr, &parms->alldone, &parms->middone);

    return HPX_SUCCESS;
  }

  /// Action to perform rank-local work
  ///
  /// This action is the target of a broadcast. Each rank will compute the
  /// DAG and start the evaluation proper. The passed in data contains a
  /// few LCOs used for synchronization.
  ///
  /// \param data - the action arguments
  /// \param msg_size - the size of the input data
  ///
  /// \returns - HPX_SUCCESS
  static int evaluate_rank_local_handler(EvaluateParams *parms,
                                         size_t msg_size) {
    // Generate the table
    RankWise<dualtree_t> global_tree{parms->rwaddr};
    auto tree = global_tree.here();
    double domain_size = tree->domain()->size();
    size_t n_params = (msg_size - sizeof(EvaluateParams)) / sizeof(double);
    std::vector<double> kernel_params(parms->kernelparams,
                                      &parms->kernelparams[n_params]);
    expansion_t::update_table(parms->n_digits, domain_size, kernel_params);
    tree->set_method(parms->method);

    // Get ready to evaluate
    // BEGIN DISTRIBUTE
    hpx_time_t distribute_begin = hpx_time_now();
    DAG *dag = tree->create_DAG();
if (hpx_get_my_rank() == 0) {
  auto sout = dag->min_max_out_degree_S();
  auto tin = dag->min_max_in_degree_T();
  auto siout = dag->min_max_out_degree_SI();
  auto tiin = dag->min_max_in_degree_TI();
  fprintf(stdout, "S: out: %d %d\n", sout.first, sout.second);
  fprintf(stdout, "SI: out: %d %d\n", siout.first, sout.second);
  fprintf(stdout, "TI: in: %d %d\n", tiin.first, tiin.second);
  fprintf(stdout, "T: in: %d %d\n", tin.first, tin.second);
}
    parms->distro.compute_distribution(*dag);
    hpx_time_t distribute_end = hpx_time_now();
    double distribute_deltat = hpx_time_diff_us(distribute_begin,
                                                distribute_end);
    //fprintf(stdout, "Evaluate: DAG creation and distribution: %lg [us]\n",
    //        distribute_deltat);
    // END DISTRIBUTE

    // BEGIN ALLOCATE
    hpx_time_t allocate_begin = hpx_time_now();
    tree->create_expansions_from_DAG(parms->rwaddr);

    // NOTE: the previous has to finish for the following. So the previous
    // is a synchronous operation. The next three, however, are not. They all
    // get their work going when they come to it and then they return.
    hpx_lco_and_set(parms->middone, HPX_NULL);
    hpx_lco_wait(parms->middone);
    hpx_time_t allocate_end = hpx_time_now();
    double allocate_deltat = hpx_time_diff_us(allocate_begin, allocate_end);
    //fprintf(stdout, "Evaluate: LCO allocation: %lg [us]\n", allocate_deltat);
    // END ALLOCATE


    // BEGIN EVALUATE
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
    tree->setup_edge_lists(dag);
    tree->start_DAG_evaluation(global_tree);
    hpx_addr_t heredone = tree->setup_termination_detection(dag);
    hpx_lco_wait(heredone);
    hpx_time_t evaluate_end = hpx_time_now();
    double evaluate_deltat = hpx_time_diff_us(evaluate_begin, evaluate_end);
    //fprintf(stdout, "Evaluate: DAG evaluation: %lg [us]\n", evaluate_deltat);

#ifdef DASHMM_INSTRUMENTATION
    libhpx_inst_phase_end();
#endif
    // END EVALUATE

    fprintf(stdout, "Evalute: %d - C/D %lg - A %lg - E %lg\n",
            hpx_get_my_rank(),
            distribute_deltat, allocate_deltat, evaluate_deltat);

    // Delete some local stuff
    hpx_lco_delete_sync(heredone);
    tree->destroy_DAG_LCOs(*dag);
    delete dag;

    // Mark that we have finished the rank-local work
    hpx_lco_and_set(parms->alldone, HPX_NULL);

    return HPX_SUCCESS;
  }

  /// Action that cleans up the tree
  ///
  /// This is called on a single locality, and will clean up the rest of the
  /// allocated resources for this evaluation. This action also exits the
  /// current HPX-5 epoch.
  ///
  /// \param rwaddr - global address of the DualTree
  /// \param alldone - global address of completion detection LCO
  /// \param middone - global address of synchronization LCO
  ///
  /// \returns HPX_SUCCESS
  static int evaluate_cleanup_handler(hpx_addr_t rwaddr, hpx_addr_t alldone,
                                      hpx_addr_t middone) {
    hpx_lco_delete_sync(alldone);
    hpx_lco_delete_sync(middone);

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
                    template <typename, typename> class> class M>
hpx_action_t Evaluator<S, T, E, M>::evaluate_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Evaluator<S, T, E, M>::evaluate_rank_local_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Evaluator<S, T, E, M>::evaluate_cleanup_ = HPX_ACTION_NULL;


} // namespace dashmm


#endif // __DASHMM_EVALUATOR_H__
