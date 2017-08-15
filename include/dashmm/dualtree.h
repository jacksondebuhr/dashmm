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


#ifndef __DASHMM_DUALTREE_H__
#define __DASHMM_DUALTREE_H__


/// \file
/// \brief DualTree related types


// C library
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cassert>

// C++ library
#include <algorithm>
#include <functional>
#include <vector>

// HPX-5
#include "hpx/hpx.h"

// DASHMM
#include "dashmm/array.h"
#include "dashmm/dag.h"
#include "dashmm/domaingeometry.h"
#include "dashmm/expansionlco.h"
#include "dashmm/hilbert.h"
#include "dashmm/index.h"
#include "dashmm/node.h"
#include "dashmm/point.h"
#include "dashmm/rankwise.h"
#include "dashmm/reductionops.h"
#include "dashmm/tree.h"


namespace dashmm {


/// Handle to a tree.
///
/// This is an opaque object that should only be used as a handle by the
/// users of DASHMM. Note that this is a handle, and should be freely
/// copied as it merely references the underlying data, and does not own that
/// data.
using DualTreeHandle = hpx_addr_t;


// Forward declare registrar so we can become friends
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class> class Method>
class DualTreeRegistrar;


/// The DualTree organizes the source and target tree and handles common work
///
/// The DualTree manages all work that instersects between the two trees.
/// This object stores the pointers to the source and target trees.
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class> class Method>
class DualTree {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, Method>;
  using targetlco_t = TargetLCO<Source, Target, Expansion, Method>;
  using method_t = Method<Source, Target, Expansion>;
  using sourcenode_t = Node<Source>;
  using targetnode_t = Node<Target>;
  using sourceref_t = ArrayRef<Source>;
  using targetref_t = ArrayRef<Target>;
  using sourcetree_t = Tree<Source, Target, Source>;
  using targettree_t = Tree<Source, Target, Target>;
  using dualtree_t = DualTree<Source, Target, Expansion, Method>;

  /// Construction is always default
  DualTree()
    : domain_{}, refinement_limit_{1}, unif_level_{1}, dim3_{8},
      unif_count_{HPX_NULL}, unif_count_value_{nullptr},
      method_{}, source_tree_{nullptr},
      target_tree_{nullptr} { }

  /// We delete the copy constructor and copy assignement operator.
  DualTree(const dualtree_t &other) = delete;
  dualtree_t &operator=(const dualtree_t &other) = delete;

  /// Return the number of source points in the uniform grid
  ///
  /// This routine returns the address of the given information,
  ///
  /// \param i - the uniform grid node in question
  int *unif_count_src(size_t i = 0) const {return &unif_count_value_[i];}

  /// Return the number of target points in the uniform grid
  ///
  /// This routine returns the address of the given information,
  ///
  /// \param i - the uniform grid node in question
  int *unif_count_tar(size_t i = 0) const {
    return &unif_count_value_[i + dim3_];
  }

  /// Return the DomainGeometry for this DualTree.
  const DomainGeometry *domain() const {return &domain_;}

  /// Return the refinement limit used to build the tree.
  int refinement_limit() const {return refinement_limit_;}

  /// Return the method this object will use for DAG operations.
  const method_t &method() const {return method_;}

  /// Set the method this object will use for DAG operations.
  void set_method(const method_t &method) {method_ = method;}

  /// Lookup target LCO address
  ///
  /// \param idx - index of LCO to look up
  /// \param op - operation of the edge leading to the LCO in question
  ///
  /// \returns - global address of LCO
  hpx_addr_t lookup_lco_addx(Index idx, Operation op) {
    bool search_source{true};
    switch(op) {
      case Operation::Nop:
        assert(0 && "Problem in address search");
        break;
      case Operation::StoM:
        break;
      case Operation::StoL:
        search_source = false;
        break;
      case Operation::MtoM:
        break;
      case Operation::MtoL:
        search_source = false;
        break;
      case Operation::LtoL:
        search_source = false;
        break;
      case Operation::MtoT:
        search_source = false;
        break;
      case Operation::LtoT:
        search_source = false;
        break;
      case Operation::StoT:
        search_source = false;
        break;
      case Operation::MtoI:
        break;
      case Operation::ItoI:
        search_source = false;
        break;
      case Operation::ItoL:
        search_source = false;
        break;
    }

    if (search_source) {
      return source_tree_.here()->lookup_lco_addx(idx, op);
    } else {
      return target_tree_.here()->lookup_lco_addx(idx, op);
    }
  }

  /// Destroy any allocated memory associated with this DualTree.
  void clear_data() {
    hpx_addr_t clear_done = hpx_lco_and_new(2);
    assert(clear_done != HPX_NULL);
    sourcetree_t::delete_tree(clear_done, source_tree_, dim3_, rank_map_);
    targettree_t::delete_tree(clear_done, target_tree_, dim3_, rank_map_);
    hpx_lco_wait(clear_done);
    hpx_lco_delete_sync(clear_done);

    // the unif counts
    delete [] unif_count_value_;
    delete [] rank_map_;
  }

  /// Return the rank owning the given unif grid node
  int rank_of_unif_grid(int idx) const {return rank_map_[idx];}

  // TODO: Get this out of DualTree
  /// Create the DAG for this tree using the method specified for this object.
  ///
  /// This will allocate and collect the DAG nodes into the returned object.
  ///
  /// \returns - the resulting DAG.
  DAG *create_DAG() {
    // Do work on the source tree
    hpx_addr_t sdone = hpx_lco_and_new(1);
    assert(sdone != HPX_NULL);

    dualtree_t *thetree = this;
    auto stree = source_tree_.here();
    // TODO: I am not happy with this. We should not have to ask for the root
    // from out here. But then what happens with apply? We do this simple thing
    // for now, and when we pull this unrelated thing out, we can fix this.
    sourcenode_t *narg = stree->root();
    hpx_call(HPX_HERE, source_apply_method_, HPX_NULL,
             &thetree, &narg, &sdone);

    hpx_lco_wait(sdone);
    hpx_lco_delete_sync(sdone);

    // Do work on the target tree
    hpx_addr_t tdone = hpx_lco_and_new(1);
    assert(tdone != HPX_NULL);

    std::vector<sourcenode_t *> *consider = new std::vector<sourcenode_t *>{};
    consider->push_back(narg);
    auto ttree = target_tree_.here();
    // TODO: I am not happy with this. We should not have to ask for the root
    // from out here. But then what happens with apply? We do this simple thing
    // for now, and when we pull this unrelated thing out, we can fix this.
    targetnode_t *trgaddx = ttree->root();
    hpx_call(HPX_HERE, target_apply_method_, HPX_NULL,
             &thetree, &trgaddx, &consider, &same_sandt_, &tdone);

    hpx_lco_wait(tdone);
    hpx_lco_delete_sync(tdone);

    DAG *retval = collect_DAG_nodes();
    return retval;
  }

  // TODO: Get this out of DualTree
  /// Traverse the tree and collect the DAG nodes
  ///
  /// This will collect the nodes of the DAG into three groups: the source
  /// nodes of the DAG, the target nodes of the DAG, and all other nodes.
  /// The source and target nodes are the input and output nodes respectively.
  /// The remainder are those nodes containing an intermediate computation.
  ///
  /// This is a synchronous operation.
  ///
  /// \param dag - a DAG object to be populated
  DAG *collect_DAG_nodes() {
    DAG *retval = new DAG{};

    auto stree = source_tree_.here();
    auto ttree = target_tree_.here();

    // TODO: same complaint with root() here...
    collect_DAG_nodes_from_S_node(stree->root(), retval->source_leaves,
                                  retval->source_nodes);
    collect_DAG_nodes_from_T_node(ttree->root(), retval->target_leaves,
                                  retval->target_nodes);

    // TODO: Note that these are non-binding requests, but this is the most
    // clear we can write this.
    retval->source_leaves.shrink_to_fit();
    retval->source_nodes.shrink_to_fit();
    retval->target_nodes.shrink_to_fit();
    retval->target_leaves.shrink_to_fit();

    return retval;
  }

  // TODO: Get this out of DualTree
  /// Create the LCOs from the DAG
  ///
  /// This will traverse the source and target tree creating any needed
  /// expansion LCOs. Further, it will create the target LCOs in the target
  /// tree.
  ///
  /// This is a synchronous operation.
  void create_expansions_from_DAG(hpx_addr_t rwtree) {
    hpx_addr_t done = hpx_lco_and_new(2);
    assert(done != HPX_NULL);

    // TODO: same complaint about root()...
    dualtree_t *argthis = this;
    sourcenode_t *srcaddx = source_tree_.here()->root();
    hpx_call(HPX_HERE, create_S_expansions_from_DAG_, HPX_NULL,
             &done, &argthis, &srcaddx, &rwtree);
    targetnode_t *trgaddx = target_tree_.here()->root();
    hpx_call(HPX_HERE, create_T_expansions_from_DAG_, HPX_NULL,
             &done, &argthis, &trgaddx, &rwtree);

    hpx_lco_wait(done);
    hpx_lco_delete_sync(done);
  }

  // TODO: Get this out of DualTree
  /// Sets up termination detection for a DASHMM evaluation
  ///
  /// This is an asynchronous operation. The returned LCO becomes the
  /// responsibility of the caller.
  ///
  /// \param targets - vector of target nodes in the DAG
  /// \param internals - vector of internal nodes in the DAG
  ///
  /// \returns - LCO that will signal that all targets have completed their
  ///            computation.
  hpx_addr_t setup_termination_detection(DAG *dag) {
    size_t n_targs = dag->target_leaves.size();
    size_t n_tinternal = dag->target_nodes.size();
    size_t n_sinternal = dag->source_nodes.size();

    hpx_addr_t retval = hpx_lco_and_new(n_targs + n_tinternal + n_sinternal);
    assert(retval != HPX_NULL);

    std::vector<DAGNode *> *argaddx = &dag->target_leaves;
    std::vector<DAGNode *> *tinternalsaddx = &dag->target_nodes;
    std::vector<DAGNode *> *sinternalsaddx = &dag->source_nodes;
    hpx_call(HPX_HERE, termination_detection_, HPX_NULL, &retval,
             &argaddx, &n_targs, &tinternalsaddx, &n_tinternal,
             &sinternalsaddx, &n_sinternal);

    return retval;
  }

  // TODO: Get this out of DualTree
  /// Initiate the DAG evaluation
  ///
  /// This starts the work of the evalution by starting the S->* work at the
  /// source nodes of the DAG.
  ///
  /// This is an asynchronous operation. The termination detection cannot
  /// possibly trigger before this is completed, so waiting on the termination
  /// of the full evaluation implicitly waits on this operation.
  ///
  /// \param global_tree - the Dual Tree
  void start_DAG_evaluation(RankWise<dualtree_t> &global_tree, DAG *dag) {
    hpx_addr_t rwaddr = global_tree.data();

    // The DAG nodes are sorted local vs not. Find the partition point
    int rank = hpx_get_my_rank();
    auto partition_point = std::partition_point(
      dag->source_leaves.begin(), dag->source_leaves.end(),
      [&rank](const DAGNode *a) -> bool {
        return a->locality == rank;
      });
    size_t count = partition_point - dag->source_leaves.begin();

    // Find chunk size
    int n_workers = hpx_get_num_threads();
    size_t delta = count / n_workers;
    if (!delta) {
      // If there are so few to spawn as to make delta < 1, then just do one
      // chunk. This should be rare.
      delta = count;
    }

    // Spawn chunks
    DAGNode **data = dag->source_leaves.data();
    size_t parc_size = sizeof(hpx_addr_t) + 2 * sizeof(DAGNode *);
    for (size_t start = 0; start < count; start += delta) {
      size_t end = start + delta;
      if (end > count) {
        end = count;
      }
      DAGNode **first = &data[start];
      DAGNode **last = &data[end];

      hpx_parcel_t *parc = hpx_parcel_acquire(nullptr, parc_size);
      hpx_parcel_set_target(parc, HPX_HERE);
      hpx_parcel_set_action(parc, instigate_dag_eval_);
      hpx_parcel_set_args(parc, &rwaddr, &first, &last);
      hpx_parcel_send_sync(parc);
    }
  }

  // TODO: Get this out of DualTree
  /// Destroys the LCOs associated with the DAG
  ///
  /// This is a synchronous operation. This destroys not only the expansion
  /// LCOs, but also the target LCOs.
  ///
  /// \param dag - the DAG
  void destroy_DAG_LCOs(DAG &dag) {
    hpx_addr_t done = hpx_lco_and_new(3);
    assert(done != HPX_NULL);

    DAGNode **data = dag.target_leaves.data();
    size_t n_data = dag.target_leaves.size();
    int type = 1;
    hpx_call(HPX_HERE, destroy_DAG_LCOs_, done, &data, &n_data, &type);

    data = dag.target_nodes.data();
    n_data = dag.target_nodes.size();
    type = 0;
    hpx_call(HPX_HERE, destroy_DAG_LCOs_, done, &data, &n_data, &type);

    data = dag.source_nodes.data();
    n_data = dag.source_nodes.size();
    type = 0;
    hpx_call(HPX_HERE, destroy_DAG_LCOs_, done, &data, &n_data, &type);

    hpx_lco_wait(done);
    hpx_lco_delete_sync(done);
  }


  /// Create the basic data for a distributed tree for use with DASHMM
  ///
  /// This will compute the domain from the given source and target points and
  /// will create the basic data for a distributed tree with that domain. The
  /// return object is a RankWise object storing each copy of the local tree.
  ///
  /// This call is synchronous, and should be called from inside an HPX thread.
  /// Further, this is to be called in a diffusive style; only a single thread
  /// should call this function.
  ///
  /// \param threshold - the partitioning threshold for the tree
  /// \param sources - the source data
  /// \param targets - the target data
  ///
  /// \returns - the RankWise object containing the dual tree
  static RankWise<dualtree_t> create(int threshold, Array<Source> sources,
                                     Array<Target> targets) {
    bool same_sandt{false};
    if (sources.data() == targets.data()) {
      same_sandt = true;
    }
    hpx_addr_t domain_geometry = compute_domain_geometry(sources, targets,
                                                         same_sandt);
    RankWise<dualtree_t> retval = setup_basic_data(threshold, domain_geometry,
                                                   same_sandt, sources,
                                                   targets);
    hpx_lco_delete_sync(domain_geometry);
    return retval;
  }

  /// Partition the tree
  ///
  /// This will do the bulk of the work for partitioning and creating the
  /// distributed tree. After the call to this routine, the source and target
  /// data will be redistributed to make evaluation easier for DASHMM. Behind
  /// the scenes, the local segments of the arrays inside sources and targets
  /// will have been replaced.
  ///
  /// This routine should be called from inside an HPX thread. Further, this is
  /// to be called in a diffusive style; this should be called from only a
  /// single thread. This routine will handle involving all ranks in the
  /// system.
  ///
  /// \param global_tree - an object previously initialized with create()
  /// \param sources - the source data
  /// \param targets - the target data
  ///
  /// \returns - an LCO indication completion of the partitioning.
  static hpx_addr_t partition(RankWise<dualtree_t> global_tree) {
    hpx_addr_t retval = hpx_lco_future_new(0);
    assert(retval != HPX_NULL);

    hpx_addr_t tree_gas = global_tree.data();
    hpx_bcast_lsync(create_dual_tree_, retval, &tree_gas);

    return retval;
  }

  /// Destroy a distributed tree.
  ///
  /// This cleans up all allocated resources used by the DualTree.
  ///
  /// This should be called from an HPX thread, in a diffusive style.
  ///
  /// \param global_tree - the distributed tree
  static void destroy(RankWise<dualtree_t> &global_tree) {
    hpx_addr_t rwtree = global_tree.data();
    hpx_bcast_rsync(finalize_partition_, &rwtree);

    auto tree = global_tree.here();
    hpx_lco_delete_sync(tree->unif_count_);

    tree->source_tree_.destroy();
    tree->target_tree_.destroy();
    global_tree.destroy();
  }


 private:
  friend class DualTreeRegistrar<Source, Target, Expansion, Method>;

  // TODO: Get the out of DualTree
  /// Edge record for DAG instigation
  struct DAGInstigationRecord {
    Operation op;
    hpx_addr_t target;
    Index idx;
  };

  /// Action to set the domain geometry given the sources and targets
  ///
  /// This action is the target of a broadcast, and computes the domain for
  /// the local sources and targets. The result is then given to a reduction
  /// LCO which reduces each rank's portion.
  ///
  /// \param sources_gas - the global address of the source data
  /// \param targets_gas - the global address of the target data
  /// \param domain_geometry - a reduction LCO to which the local domain is sent
  /// \param same_sandt - is this an S==T DualTree
  ///
  /// \returns HPX_SUCCESS
  static int set_domain_geometry_handler(hpx_addr_t sources_gas,
                                         hpx_addr_t targets_gas,
                                         hpx_addr_t domain_geometry,
                                         int same_sandt) {
    Array<source_t> sources{sources_gas};
    sourceref_t src_ref = sources.ref();
    Source *s = src_ref.data();

    double var[6] = {1e50, -1e50, 1e50, -1e50, 1e50, -1e50};

    for (size_t i = 0; i < src_ref.n(); ++i) {
      var[0] = fmin(var[0], s[i].position.x());
      var[1] = fmax(var[1], s[i].position.x());
      var[2] = fmin(var[2], s[i].position.y());
      var[3] = fmax(var[3], s[i].position.y());
      var[4] = fmin(var[4], s[i].position.z());
      var[5] = fmax(var[5], s[i].position.z());
    }

    // Only do the targets if they are different from the sources
    if (!same_sandt) {
      Array<target_t> targets{targets_gas};
      targetref_t trg_ref = targets.ref();
      Target *t = trg_ref.data();

      for (size_t i = 0; i < trg_ref.n(); ++i) {
        var[0] = fmin(var[0], t[i].position.x());
        var[1] = fmax(var[1], t[i].position.x());
        var[2] = fmin(var[2], t[i].position.y());
        var[3] = fmax(var[3], t[i].position.y());
        var[4] = fmin(var[4], t[i].position.z());
        var[5] = fmax(var[5], t[i].position.z());
      }
    }

    hpx_lco_set_lsync(domain_geometry, sizeof(double) * 6, var, HPX_NULL);

    return HPX_SUCCESS;
  }

  /// Operation implementing the identity for the domain reduction
  ///
  /// \param values - the LCO's data
  /// \param UNUSED - the size of the data
  static void domain_geometry_init_handler(double *values,
                                           const size_t UNUSED) {
    values[0] = 1e50; // xmin
    values[1] = -1e50; // xmax
    values[2] = 1e50; // ymin
    values[3] = -1e50; // ymax
    values[4] = 1e50; // zmin
    values[5] = -1e50; // zmax
  }

  /// Operation implementing the reduction for the domain reduction
  ///
  /// \param lhs - the LCO data
  /// \param rhs - the set data
  /// \param UNUSED - the size of the set data
  static void domain_geometry_op_handler(double *lhs,
                                         double *rhs,
                                         size_t UNUSED) {
    lhs[0] = fmin(lhs[0], rhs[0]);
    lhs[1] = fmax(lhs[1], rhs[1]);
    lhs[2] = fmin(lhs[2], rhs[2]);
    lhs[3] = fmax(lhs[3], rhs[3]);
    lhs[4] = fmin(lhs[4], rhs[4]);
    lhs[5] = fmax(lhs[5], rhs[5]);
  }

  /// Compute the bounding box for the given source and target points.
  ///
  /// This will return an LCO which will contain six doubles, (xmin, xmax,
  /// ymin, ymax, zmin, zmax) This is an asynchronous call, it will return
  /// as soon as the work is scheduled.
  ///
  /// \param sources - the source array
  /// \param targets - the target array
  /// \param same_sandt - is S == T for this tree
  ///
  /// \returns - address of an LCO containing the reduced domain
  static hpx_addr_t compute_domain_geometry(Array<Source> sources,
                                            Array<Target> targets,
                                            bool same_sandt) {
    // Create a reduction LCO
    hpx_addr_t domain_geometry =
      hpx_lco_reduce_new(hpx_get_num_ranks(), sizeof(double) * 6,
                         domain_geometry_init_,
                         domain_geometry_op_);

    // Launch the reduction actions
    hpx_addr_t sglob = sources.data();
    hpx_addr_t tglob = targets.data();
    int ssat = (same_sandt ? 1 : 0);
    hpx_bcast_lsync(set_domain_geometry_, HPX_NULL,
                    &sglob, &tglob, &domain_geometry, &ssat);

    return domain_geometry;
  }

  /// Action to perform initializtion of basic data for the local tree
  ///
  /// This is the target of a broadcast, and it sets various data about the
  /// tree.
  ///
  /// \param rwdata - the global address of the global tree
  /// \param count - an LCO in which the uniform grid counting is reduced
  /// \param limit - the partitioning threshold for the tree
  /// \param domain_geometry - the LCO in which the domain is reduced
  /// \param same_sandt - is S == T for this tree
  /// \param source_gas - the source records
  /// \param target_gas - the target records
  /// \param stree - the global address of the source tree
  /// \param ttree - the global address of the target tree
  ///
  /// \returns - HPX_SUCCESS
  static int init_partition_handler(hpx_addr_t rwdata,
                                    hpx_addr_t count,
                                    int limit,
                                    hpx_addr_t domain_geometry,
                                    int same_sandt,
                                    hpx_addr_t source_gas,
                                    hpx_addr_t target_gas,
                                    hpx_addr_t stree,
                                    hpx_addr_t ttree) {
    RankWise<dualtree_t> global_tree{rwdata};
    auto tree = global_tree.here();

    int num_ranks = hpx_get_num_ranks();
    tree->unif_level_ = ceil(log(num_ranks) / log(8)) + 1;
    tree->dim3_ = pow(8, tree->unif_level_);
    tree->unif_count_ = count;
    tree->refinement_limit_ = limit;
    tree->source_tree_ = RankWise<sourcetree_t>{stree};
    tree->target_tree_ = RankWise<targettree_t>{ttree};
    tree->same_sandt_ = same_sandt;
    tree->source_gas = source_gas;
    tree->target_gas = target_gas;

    // Call out to tree setup stuff
    hpx_addr_t setup_done = hpx_lco_and_new(2);
    assert(setup_done != HPX_NULL);

    {
      auto source_tree = tree->source_tree_.here();
      *source_tree = sourcetree_t{};
      source_tree->setupBasics(setup_done, tree->unif_level_);
    }

    {
      auto target_tree = tree->target_tree_.here();
      *target_tree = targettree_t{};
      target_tree->setupBasics(setup_done, tree->unif_level_);
    }

    // We here allocate space for the result of the counting
    tree->unif_count_value_ = new int[tree->dim3_ * 2]();

    // Setup domain_
    double var[6];
    hpx_lco_get(domain_geometry, sizeof(double) * 6, &var);
    double length = fmax(var[1] - var[0],
                         fmax(var[3] - var[2], var[5] - var[4]));
    DomainGeometry geo{Point{(var[1] + var[0] - length) / 2,
                             (var[3] + var[2] - length) / 2,
                             (var[5] + var[4] - length) / 2}, length};
    tree->domain_ = geo;

    // Wait for setup to be done
    hpx_lco_wait(setup_done);
    hpx_lco_delete_sync(setup_done);

    return HPX_SUCCESS;
  }

  /// Allocate and setup a Dual Tree
  ///
  /// This will both allocate and setup a dual tree.
  ///
  /// \param threshold - the partitioning threshold
  /// \param domain_geometry - an LCO into which the domain is reduced
  /// \param same_sandt - is S == T for this tree
  /// \param sources - the source Array
  /// \param targets - the target Array
  ///
  /// \returns - the Dual Tree
  static RankWise<dualtree_t> setup_basic_data(int threshold,
                                               hpx_addr_t domain_geometry,
                                               bool same_sandt,
                                               Array<source_t> sources,
                                               Array<target_t> targets) {
    RankWise<dualtree_t> retval{};
    retval.allocate();
    assert(retval.valid());

    RankWise<sourcetree_t> stree{};
    stree.allocate();
    RankWise<targettree_t> ttree{};
    ttree.allocate();
    assert(stree.valid() && ttree.valid());

    // Now the single things are created.
    int num_ranks = hpx_get_num_ranks();
    int level = ceil(log(num_ranks) / log(8)) + 1;
    int dim3 = pow(8, level);
    hpx_addr_t ucount = hpx_lco_reduce_new(num_ranks, sizeof(int) * (dim3 * 2),
                                           int_sum_ident_op,
                                           int_sum_op);
    hpx_addr_t rwdata = retval.data();
    int ssat = (same_sandt ? 1 : 0);
    hpx_addr_t sgas = sources.data();
    hpx_addr_t tgas = targets.data();
    hpx_addr_t stree_addx = stree.data();
    hpx_addr_t ttree_addx = ttree.data();
    hpx_bcast_rsync(init_partition_, &rwdata, &ucount, &threshold,
                    &domain_geometry, &ssat, &sgas, &tgas, &stree_addx,
                    &ttree_addx);

    return retval;
  }

  /// Distribute the uniforn level nodes
  ///
  /// Given the number of sources and targets in each node at the uniform level
  /// of the tree, this will distribute those uniform level nodes amoung the
  /// given number of ranks. This will allocate an array and return that
  /// array to the caller; the caller assumes ownership of that array.
  ///
  /// The input counts are provided as an array of integers that are back to
  /// back, and have a length @p len, with the source counts per node in the
  /// first half of @p global and the target counts per node in the second
  /// half of @p global.
  ///
  /// The mapping of the nodes to index in the provided @p global and the
  /// returned array is according to a simple index. At the uniform level, if
  /// the uniform level is l, there are N = 2^l nodes along each length of the
  /// cubical domain. So indexed from the low corner, these nodes can be
  /// described with three indices (ix, iy, iz). The order of each node in the
  /// input and output data is computed with: ix + iy * N  + iz * N * N.
  ///
  /// \param num_ranks - the number of ranks over which to distribute the nodes
  /// \param global - the source and target counts per node
  /// \param len - the number of uniform level nodes.
  /// \param lvl - the uniform level.
  ///
  /// \returns - an array assigning each node to a rank. Caller assumes
  ///            ownership of this array, and should destroy it when done with
  ///            the data.
  int *distribute_points(int num_ranks, const int *global, int len, int lvl) {
    return distribute_points_hilbert(num_ranks, global, len, lvl);
  }

  /// Count and sort the local points
  ///
  /// This will assign the local points to the uniform grid, and it will also
  /// rearrange them according to which node of the uniform grid. This will
  /// ultimately return the local counts per uniform grid node which will later
  /// be combined into a global count. The returned counts are allocated in
  /// this routine; the caller assumes ownership of the returned array.
  ///
  /// \param tree - the dual tree
  /// \param p_s - the source data
  /// \param n_sources - the number of sources
  /// \param p_t - the target data
  /// \param n_targets - the number of targets
  /// \param local_offset_s - array that will be allocated and filled with the
  ///                         local offsets for the sources
  /// \param local_offset_t - array that will be allocated and filled with the
  ///                         local offsets for the targets
  ///
  /// \returns - the local counts per uniform grid node
  static int *sort_local_points(DualTree *tree,
                                source_t *p_s,
                                int n_sources,
                                target_t *p_t,
                                int n_targets,
                                int **local_offset_s,
                                int **local_offset_t) {
    int *local_count = new int[tree->dim3_ * 2]();
    int *local_scount = local_count;
    int *local_tcount = &local_count[tree->dim3_];

    int *gid_of_sources = new int[n_sources]();
    int *gid_of_targets = new int[n_targets]();

    // Assign points to the grid
    hpx_addr_t assign_done = hpx_lco_and_new(2);
    assert(assign_done != HPX_NULL);
    sourcetree_t::assign_points_to_unif_grid(assign_done, p_s, n_sources,
                                             &tree->domain_, tree->unif_level_,
                                             gid_of_sources, local_scount);
    targettree_t::assign_points_to_unif_grid(assign_done, p_t, n_targets,
                                             &tree->domain_, tree->unif_level_,
                                             gid_of_targets, local_tcount);
    hpx_lco_wait(assign_done);
    hpx_lco_delete_sync(assign_done);

    // Exchange counts
    hpx_lco_set(tree->unif_count_, sizeof(int) * tree->dim3_ * 2, local_count,
                HPX_NULL, HPX_NULL);

    // Group points on the same grid
    hpx_addr_t group_done = hpx_lco_and_new(2);
    assert(group_done != HPX_NULL);
    sourcetree_t::group_points_on_unif_grid(group_done, p_s, n_sources,
                                            tree->dim3_, gid_of_sources,
                                            local_scount, local_offset_s);
    if (!tree->same_sandt_) {
      targettree_t::group_points_on_unif_grid(group_done, p_t, n_targets,
                                              tree->dim3_, gid_of_targets,
                                              local_tcount, local_offset_t);
    } else {
      hpx_lco_and_set(group_done, HPX_NULL);
    }
    hpx_lco_wait(group_done);
    hpx_lco_delete_sync(group_done);

    delete [] gid_of_sources;
    delete [] gid_of_targets;

    return local_count;
  }

  /// Merge incoming points into the sorted list and spawns adapative partition
  ///
  /// This action is the 'far side' of the send points message. This will
  /// read the incoming points and merge them with the local data. If the
  /// uniform grid node to which they belong has received all of the points it
  /// is waiting for, this routine will then spawn the adaptive partitioning
  /// of that branch.
  ///
  /// This is a marshalled action, and so the message data is rather opaque.
  ///
  /// \param args - a buffer containing the incoming message.
  /// \parma UNUSED - the size of the message.
  static int recv_points_handler(void *args, size_t UNUSED) {
    hpx_addr_t *rwarg = static_cast<hpx_addr_t *>(args);
    RankWise<dualtree_t> global_tree{*rwarg};
    auto local_tree = global_tree.here();

    Serializer *source_manager{nullptr};
    Serializer *target_manager{nullptr};
    {
      Array<source_t> sarr{rwarg[1]};
      source_manager = sarr.get_manager();
      Array<target_t> tarr{rwarg[2]};
      target_manager = tarr.get_manager();
    }

    // Wait until the buffer is allocated before merging incoming messages
    // We could do this as a call when, but then we need to be aware of the
    // addresses for every rank's LCO. For now, we do this, as it is simpler.
    auto stree = local_tree->source_tree_.here();
    auto ttree = local_tree->target_tree_.here();
    stree->waitForUnifDone();
    ttree->waitForUnifDone();

    int *meta = reinterpret_cast<int *>(static_cast<char *>(args)
                                        + 3 * sizeof(hpx_addr_t));
    int range = meta[0];
    int recv_ns = meta[1];
    int recv_nt = meta[2];
    int *count_s = &meta[3]; // Used only if recv_ns > 0
    int *count_t = count_s + range * (recv_ns > 0); // Used only if recv_nt > 0
    char *recv_s = static_cast<char *>(args) + sizeof(int) * 3 +
                                3 * sizeof(hpx_addr_t) +
                                sizeof(int) * range * (recv_ns > 0) +
                                sizeof(int) * range * (recv_nt > 0);
    // Note that there are two copies now, the first is as the parcel data is
    // deserialized, and the second is to add it to the buffer. This is
    // essentially unavoidable.
    source_t *sources = new source_t[recv_ns];
    for (int i = 0; i < recv_ns; ++i) {
      recv_s = (char *)source_manager->deserialize(recv_s, &sources[i]);
    }
    target_t *targets{nullptr};
    if (recv_nt) {
      targets = new target_t[recv_nt];
      for (int i = 0; i < recv_nt; ++i) {
        recv_s = (char *)target_manager->deserialize(recv_s, &targets[i]);
      }
    }

    hpx_addr_t done = hpx_lco_and_new(range * 2);

    int myrank = hpx_get_my_rank();
    if (recv_ns) {
      int offset{0};
      int current_box{0};
      for (int i = 0; i < local_tree->dim3_; ++i) {
        if (local_tree->rank_of_unif_grid(i) == myrank) {
          int incoming_ns = count_s[current_box];
          if (incoming_ns) {
            stree->mergePoints(done, &sources[offset], i, incoming_ns,
                               local_tree->domain(),
                               local_tree->refinement_limit());
            offset += incoming_ns;
          } else {
            hpx_lco_and_set(done, HPX_NULL);
          }
          ++current_box;
        }
      }
      // After that we should have treated the number of boxes the sender
      // computed as needing data.
      assert(current_box == range);
    } else {
      if (range) {
        hpx_lco_and_set_num(done, range, HPX_NULL);
      }
    }

    if (recv_nt) {
      int offset{0};
      int current_box{0};
      for (int i = 0; i < local_tree->dim3_; ++i) {
        if (local_tree->rank_of_unif_grid(i) == myrank) {
          int incoming_nt = count_t[current_box];
          if (incoming_nt) {
            ttree->mergePoints(done, &targets[offset], i, incoming_nt,
                               local_tree->domain(),
                               local_tree->refinement_limit());
            offset += incoming_nt;
          } else {
            hpx_lco_and_set(done, HPX_NULL);
          }
          ++current_box;
        }
      }
      // After that we should have treated the number of boxes the sender
      // computed as needing data. Further, we do this test in both places
      // in case it should be that there were no sources sent, but only
      // tagets.
      assert(current_box == range);
    } else {
      if (local_tree->same_sandt_ && recv_ns) {
        // S == T means do a special version of merge.
        int current_box{0};
        for (int i = 0; i < local_tree->dim3_; ++i) {
          if (local_tree->rank_of_unif_grid(i) == myrank) {
            int incoming_nt = count_s[current_box];
            if (incoming_nt) {
              ttree->mergePointsSameSAndT(done, i, incoming_nt, stree.local(),
                                          local_tree->domain(),
                                          local_tree->refinement_limit());
            } else {
              hpx_lco_and_set(done, HPX_NULL);
            }
            ++current_box;
          }
        }
        // Same comment.
        assert(current_box == range);
      } else {
        if (range) {
          hpx_lco_and_set_num(done, range, HPX_NULL);
        }
      }
    }

    // Wait until the data has been merged before releasing the parcel
    // NOTE: This 'done' will trigger once points are merged. The action that
    // triggers this will spawn more work, but not in a synchronous way.
    hpx_lco_wait(done);
    hpx_lco_delete_sync(done);

    delete [] sources;
    if (recv_nt) {
      delete [] targets;
    }

    return HPX_SUCCESS;
  }

  /// Send the points to the remote that will assume ownership
  ///
  /// This action sends points to remote localities that have been assigned
  /// the points during distribution.
  ///
  /// \param rank - the rank to which we are sending
  /// \param count_s - the source counts
  /// \param count_t - the target counts
  /// \param offset_s - the source offsets
  /// \param offset_t - the target offsets
  /// \param sources - the source data
  /// \param targets - the target data
  /// \param smeta - the source meta data global address
  /// \param tmeta - the target meta data global address
  /// \param rwaddr - the global address of the dual tree
  ///
  /// \returns - HPX_SUCCESS
  static int send_points_handler(int rank,
                                 int *count_s,
                                 int *count_t,
                                 int *offset_s,
                                 int *offset_t,
                                 source_t *sources,
                                 target_t *targets,
                                 hpx_addr_t smeta,
                                 hpx_addr_t tmeta,
                                 hpx_addr_t rwaddr) {
    RankWise<dualtree_t> global_tree{rwaddr};
    auto local_tree = global_tree.here();

    // Note: all the pointers are local to the calling rank.
    int range{0};
    int send_ns = 0;
    int send_nt = 0;
    for (int i = 0; i < local_tree->dim3_; ++i) {
      if (local_tree->rank_of_unif_grid(i) == rank) {
        send_ns += count_s[i];
        send_nt += count_t[i];
        ++range;
      }
    }

    // Clear out the target sends if S == T
    if (local_tree->same_sandt_) {
      send_nt = 0;
    }

    // get access to the serializers
    Serializer *source_manager{nullptr};
    Serializer *target_manager{nullptr};
    {
      Array<source_t> sarr{smeta};
      source_manager = sarr.get_manager();
      Array<target_t> tarr{tmeta};
      target_manager = tarr.get_manager();
    }

    // Parcel message length
    size_t bytes = sizeof(hpx_addr_t) * 3 + sizeof(int) * 3;
    size_t bytes_source{0};
    if (send_ns) {
      bytes += sizeof(int) * range;
      for (int i = 0; i < local_tree->dim3_; ++i) {
        if (local_tree->rank_of_unif_grid(i) == rank) {
          for (int j = 0; j < count_s[i]; ++j) {
            bytes_source += source_manager->size(&sources[offset_s[i] + j]);
          }
          bytes += bytes_source;
        }
      }
    }
    if (send_nt) {
      bytes += sizeof(int) * range;
      for (int i = 0; i < local_tree->dim3_; ++i) {
        if (local_tree->rank_of_unif_grid(i) == rank) {
          for (int j = 0; j < count_t[i]; ++j) {
            bytes += target_manager->size(&targets[offset_t[i] + j]);
          }
        }
      }
    }

    // Acquire parcel
    hpx_parcel_t *p = hpx_parcel_acquire(nullptr, bytes);
    void *data = hpx_parcel_get_data(p);
    hpx_addr_t *rwarg = static_cast<hpx_addr_t *>(data);
    rwarg[0] = rwaddr;
    rwarg[1] = smeta;
    rwarg[2] = tmeta;
    int *meta = reinterpret_cast<int *>(static_cast<char *>(data)
                                          + 3 * sizeof(hpx_addr_t));
    meta[0] = range;
    meta[1] = send_ns;
    meta[2] = send_nt;

    int *count = &meta[3];
    if (send_ns) {
      for (int i = 0; i < local_tree->dim3_; ++i) {
        if (local_tree->rank_of_unif_grid(i) == rank) {
          *count = count_s[i];
          count += 1;
        }
      }
    }
    if (send_nt) {
      for (int i = 0; i < local_tree->dim3_; ++i) {
        if (local_tree->rank_of_unif_grid(i) == rank) {
          *count = count_t[i];
          count += 1;
        }
      }
    }

    char *meta_s = static_cast<char *>(data) + 3 * sizeof(hpx_addr_t) +
      sizeof(int) * (3 + range * (send_ns > 0) + range * (send_nt > 0));
    if (send_ns) {
      for (int i = 0; i < local_tree->dim3_; ++i) {
        if (local_tree->rank_of_unif_grid(i) == rank) {
          for (int j = 0; j < count_s[i]; ++j) {
            meta_s = (char *)source_manager->serialize(
                                &sources[offset_s[i] + j],
                                meta_s);
          }
        }
      }
    }

    if (send_nt) {
      for (int i = 0; i < local_tree->dim3_; ++i) {
        if (local_tree->rank_of_unif_grid(i) == rank) {
          for (int j = 0; j < count_t[i]; ++j) {
            meta_s = (char *)target_manager->serialize(
                                &targets[offset_t[i] + j],
                                meta_s);
          }
        }
      }
    }

    hpx_parcel_set_target(p, HPX_THERE(rank));
    hpx_parcel_set_action(p, recv_points_);
    hpx_parcel_send(p, HPX_NULL);

    return HPX_SUCCESS;
  }

  /// This will prune links in the top part of the tree that are not needed
  ///
  /// After exchanging counts, there could be some nodes at the uniform level
  /// that will not have any particles under them. In this case, the links to
  /// these nodes from higher levels should be culled to avoid sending Methods
  /// into pointless work, and potentially scheduling pointless work.
  ///
  /// TODO: At the moment, this is done in serial. It might be useful to
  /// get some parallelism here. Probably the easiest is to do this work
  /// as its own task.
  void prune_topnodes() {
    int *s_counts = unif_count_src();
    auto stree = source_tree_.here();
    stree->pruneTopnodes(s_counts, dim3_, unif_level_);

    int *t_counts = unif_count_tar();
    auto ttree = target_tree_.here();
    ttree->pruneTopnodes(t_counts, dim3_, unif_level_);
  }

  /// The action responsible for dual tree partitioning
  ///
  /// This action is the target of a broadcast and manages all the work
  /// required to build the distributed trees.
  ///
  /// \param rwtree - the global address of the dual tree
  ///
  /// \return HPX_SUCCESS
  static int create_dual_tree_handler(hpx_addr_t rwtree) {
    int rank = hpx_get_my_rank();
    int num_ranks = hpx_get_num_ranks();

    RankWise<dualtree_t> global_tree{rwtree};
    auto tree = global_tree.here();
    auto sources_gas = tree->source_gas;
    auto targets_gas = tree->target_gas;

    Array<source_t> sources{sources_gas};
    Array<target_t> targets{targets_gas};

    auto stree = tree->source_tree_.here();
    auto ttree = tree->target_tree_.here();

    {
      sourceref_t src_ref = sources.ref();
      source_t *p_s = src_ref.data();
      int n_sources = src_ref.n();

      targetref_t trg_ref = targets.ref();
      target_t *p_t = trg_ref.data();
      int n_targets = trg_ref.n();

      // Assign points to uniform grid
      int *local_offset_s{nullptr};
      int *local_offset_t{nullptr};
      int *local_count = sort_local_points(&*tree, p_s, n_sources,
                                           p_t, n_targets, &local_offset_s,
                                           &local_offset_t);
      int *local_scount = local_count;
      int *local_tcount = &local_count[tree->dim3_];

      // Compute point distribution
      hpx_lco_get(tree->unif_count_, sizeof(int) * (tree->dim3_ * 2),
                  tree->unif_count_src());
      tree->rank_map_ = tree->distribute_points(num_ranks,
                                                tree->unif_count_src(),
                                                tree->dim3_,
                                                tree->unif_level_);

      // Exchange points
      stree->initPointExchange(tree->dim3_, tree->rank_map_,
                               tree->unif_count_src(), local_scount,
                               local_offset_s, tree->domain_,
                               tree->refinement_limit_, p_s);

      if (!tree->same_sandt_) {
        ttree->initPointExchange(tree->dim3_, tree->rank_map_,
                                 tree->unif_count_tar(), local_tcount,
                                 local_offset_t, tree->domain_,
                                 tree->refinement_limit_, p_t);
      } else {
        ttree->initPointExchangeSameSAndT(tree->dim3_, tree->rank_map_,
                                          tree->unif_count_tar(), local_tcount,
                                          local_offset_t, tree->domain(),
                                          tree->refinement_limit_, p_t,
                                          stree.local());
      }

      // So this one is pretty simple. It sends those points from this rank
      // going to the other rank in a parcel.
      for (int r = 0; r < num_ranks; ++r) {
        if (r != rank) {
          hpx_call(HPX_HERE, send_points_, HPX_NULL, &r, &local_scount,
                   &local_tcount, &local_offset_s, &local_offset_t,
                   &p_s, &p_t, &sources_gas, &targets_gas, &rwtree);
        }
      }

      hpx_addr_t dual_tree_complete = hpx_lco_and_new(2 * tree->dim3_);
      assert(dual_tree_complete != HPX_NULL);

      for (int i = 0; i < tree->dim3_; ++i) {
        if (tree->rank_of_unif_grid(i) == rank) {
          if (*(tree->unif_count_src(i)) == 0) {
            hpx_lco_and_set(dual_tree_complete, HPX_NULL);
          } else {
            stree->sendNode(dual_tree_complete, i, tree->source_tree_.data());
          }

          if (*(tree->unif_count_tar(i)) == 0) {
            hpx_lco_and_set(dual_tree_complete, HPX_NULL);
          } else {
            ttree->sendNode(dual_tree_complete, i, tree->target_tree_.data());
          }
        } else {
          if (*(tree->unif_count_src(i)) == 0) {
            hpx_lco_and_set(dual_tree_complete, HPX_NULL);
          } else {
            auto complco = stree->unifNodeCompletion(i);
            hpx_call_when_with_continuation(complco,
                dual_tree_complete, hpx_lco_set_action,
                complco, hpx_lco_delete_action,
                nullptr, 0);
          }

          if (*(tree->unif_count_tar(i)) == 0) {
            hpx_lco_and_set(dual_tree_complete, HPX_NULL);
          } else {
            auto complco = ttree->unifNodeCompletion(i);
            hpx_call_when_with_continuation(complco,
                dual_tree_complete, hpx_lco_set_action,
                complco, hpx_lco_delete_action,
                nullptr, 0);
          }
        }
      }

      hpx_lco_wait(dual_tree_complete);
      hpx_lco_delete_sync(dual_tree_complete);

      // This will prune pointless nodes from the top of the tree.
      tree->prune_topnodes();

      delete [] local_count;
      delete [] local_offset_s;
      delete [] local_offset_t;
    }

    // Replace segment in the array
    source_t *old_src_data = sources.replace(stree->sorted());
    if (old_src_data != nullptr) {
      delete [] old_src_data;
    }
    if (!tree->same_sandt_) {
      target_t *old_tar_data = targets.replace(ttree->sorted());
      if (old_tar_data != nullptr) {
        delete [] old_tar_data;
      }
    }

    return HPX_SUCCESS;
  }

  /// Action to destroy the tree
  ///
  /// This action is the target of a broadcast that is used to destroy the
  /// tree.
  ///
  /// \param rwtree - global address of the dual tree
  ///
  /// \returns - HPX_SUCCESS
  static int finalize_partition_handler(hpx_addr_t rwtree) {
    RankWise<dualtree_t> global_tree{rwtree};
    auto tree = global_tree.here();
    tree->clear_data();
    return HPX_SUCCESS;
  }

  // TODO: Get this out of DualTree
  /// Collect DAG nodes from source Tree Nodes
  ///
  /// \param root - tree node
  /// \param sources [out] - DAG nodes that are sources
  /// \param internals [out] - other source tree DAG nodes
  void collect_DAG_nodes_from_S_node(sourcenode_t *root,
                                     std::vector<DAGNode *> &sources,
                                     std::vector<DAGNode *> &internals) {
    for (int i = 0; i < 8; ++i) {
      if (root->child[i]) {
        collect_DAG_nodes_from_S_node(root->child[i], sources, internals);
      }
    }
    root->dag.collect_DAG_nodes(sources, internals);
  }

  // TODO: Get this out of DualTree
  /// Collect DAG nodes from target Tree Nodes
  ///
  /// \param root - tree node
  /// \param targets [out] - DAG nodes that are targets
  /// \param internals [out] - other target tree DAG nodes
  void collect_DAG_nodes_from_T_node(targetnode_t *root,
                                     std::vector<DAGNode *> &targets,
                                     std::vector<DAGNode *> &internals) {
    for (int i = 0; i < 8; ++i) {
      if (root->child[i]) {
        collect_DAG_nodes_from_T_node(root->child[i], targets, internals);
      }
    }
    root->dag.collect_DAG_nodes(targets, internals);
  }

  /// TODO: Get this out of DualTree
  /// Action to create LCOs from the DAG
  ///
  /// This will walk through the source tree and create the Expansion LCOs
  /// needed by @p node. This will spawns recursively as HPX-5 actions the
  /// full tree.
  ///
  /// \param done - LCO address for indicating completion of LCO creation
  /// \param tree - the DualTree
  /// \param node - the node under examination
  /// \param rwtree - the global address of the DualTree
  ///
  /// \returns - HPX_SUCCESS
  static int create_S_expansions_from_DAG_handler(hpx_addr_t done,
                                                  dualtree_t *tree,
                                                  sourcenode_t *node,
                                                  hpx_addr_t rwtree) {
    Point n_center = tree->domain_.center_from_index(node->idx);

    int myrank = hpx_get_my_rank();

    // create the normal expansion if needed
    if (node->dag.has_normal() && node->dag.normal()->locality == myrank) {
      expansionlco_t expand(node->dag.normal(),
                            node->idx, kSourcePrimary,
                            expansion_t::compute_scale(node->idx),
                            n_center,
                            rwtree);
      node->dag.set_normal_expansion(expand.lco());
    }

    // If there is to be an intermediate expansion, create that
    if (node->dag.has_interm() && node->dag.interm()->locality == myrank) {
      expansionlco_t intexp_lco(node->dag.interm(),
                                node->idx,
                                kSourceIntermediate,
                                expansion_t::compute_scale(node->idx),
                                n_center,
                                rwtree);
      node->dag.set_interm_expansion(intexp_lco.lco());
    }

    // spawn work at children
    int n_children{node->n_children()};

    if (n_children) {
      hpx_addr_t cdone = hpx_lco_and_new(n_children);
      assert(cdone != HPX_NULL);

      for (int i = 0; i < 8; ++i) {
        if (node->child[i] != nullptr) {
          sourcenode_t *srcaddx = node->child[i];
          hpx_call(HPX_HERE, create_S_expansions_from_DAG_, HPX_NULL,
                   &cdone, &tree, &srcaddx, &rwtree);
        }
      }

      // This will set the parent's LCO as well as delete cdone
      hpx_call_when(cdone, cdone, hpx_lco_delete_action,
                    done, nullptr, 0);
    } else {
      hpx_lco_set(done, 0, nullptr, HPX_NULL, HPX_NULL);
    }

    return HPX_SUCCESS;
  }

  // TODO: Get this out of DualTree
  /// Action to create LCOs from the DAG
  ///
  /// This will walk through the target tree and create the Expansion LCOs
  /// needed by @p node. This will spawns recursively as HPX-5 actions the
  /// full tree.
  ///
  /// \param done - LCO address for indicating completion of LCO creation
  /// \param tree - the DualTree
  /// \param node - the node under examination
  /// \param rwtree - the global address of the DualTree
  ///
  /// \returns - HPX_SUCCESS
  static int create_T_expansions_from_DAG_handler(hpx_addr_t done,
                                                  dualtree_t *tree,
                                                  targetnode_t *node,
                                                  hpx_addr_t rwtree) {
    Point n_center = tree->domain_.center_from_index(node->idx);

    int myrank = hpx_get_my_rank();

    // create the normal expansion if needed
    if (node->dag.has_normal() && node->dag.normal()->locality == myrank) {
      expansionlco_t expand(node->dag.normal(),
                            node->idx,
                            kTargetPrimary,
                            expansion_t::compute_scale(node->idx),
                            n_center,
                            rwtree);
      node->dag.set_normal_expansion(expand.lco());
    }

    // If there is to be an intermediate expansion, create that
    if (node->dag.has_interm() && node->dag.interm()->locality == myrank) {
      expansionlco_t intexp_lco(node->dag.interm(),
                                node->idx,
                                kTargetIntermediate,
                                expansion_t::compute_scale(node->idx),
                                n_center,
                                rwtree);
      node->dag.set_interm_expansion(intexp_lco.lco());
    }

    // NOTE: this spawn through the tree does not end when the tree ends.
    // Instead, we have to check if this node has a parts node in the DAG.
    // If so, this branch is done, and we need not spawn more.

    // Here is where we make the target lco if needed
    if (node->dag.has_parts()) {
      if (node->dag.parts()->locality == myrank) {
        targetlco_t tlco{node->dag.parts()->in_count(), node->parts};
        node->dag.set_targetlco(tlco.lco());
      }

      hpx_lco_set(done, 0, nullptr, HPX_NULL, HPX_NULL);
    } else {
      hpx_addr_t cdone = hpx_lco_and_new(node->n_children());
      assert(cdone != HPX_NULL);

      for (int i = 0; i < 8; ++i) {
        if (node->child[i] != nullptr) {
          targetnode_t *trgaddx = node->child[i];
          hpx_call(HPX_HERE, create_T_expansions_from_DAG_, HPX_NULL,
                   &cdone, &tree, &trgaddx, &rwtree);
        }
      }

      // This will set the parent's LCO as well as delete cdone
      hpx_call_when(cdone, cdone, hpx_lco_delete_action,
                    done, nullptr, 0);
    }

    return HPX_SUCCESS;
  }

  // TODO: Get this out of DualTree
  /// Action to start DAG evaluation work
  ///
  /// This is a recursove spawn through the source tree to begin the work of
  /// the DAG evaluation. This will cause the source DAG nodes to begin their
  /// out edges.
  ///
  /// \param rwtree - the global address of the DualTree
  /// \param first - the first DAGNode in consideration
  /// \param last - the last DAGNode in consideration
  ///
  /// \returns - HPX_SUCCESS
  static int instigate_dag_eval_handler(hpx_addr_t rwtree,
                                        DAGNode **first,
                                        DAGNode **last) {
    RankWise<dualtree_t> global_tree{rwtree};
    auto tree = global_tree.here();
    Serializer *manager{nullptr};
    {
      Array<source_t> sarr{tree->source_gas};
      manager = sarr.get_manager();
    }

    for (DAGNode **iter = first; iter != last; ++iter) {
      DAGNode *parts = *iter;
      sourcenode_t *node = static_cast<sourcenode_t *>(parts->tree_node());
      sourceref_t sources = node->parts;

      // We first sort the out edges by locality
      std::sort(parts->out_edges.begin(), parts->out_edges.end(),
                DAG::compare_edge_locality);

      // Make scratch space for the sends
      size_t source_size{0};
      auto sref = sources.data();
      for (size_t i = 0; i < sources.n(); ++i) {
        source_size += manager->size(&sref[i]);
      }
      size_t header_size = source_size + sizeof(size_t)
          + sizeof(hpx_addr_t);
      size_t total_size = header_size + sizeof(size_t)
          + parts->out_edges.size() * sizeof(DAGInstigationRecord);
      char *scratch = new char [total_size];

      // Copy source data
      {
        hpx_addr_t *scratch_rw = reinterpret_cast<hpx_addr_t *>(scratch);
        *scratch_rw = rwtree;

        size_t *scratch_n = reinterpret_cast<size_t *>(scratch
                                                       + sizeof(hpx_addr_t));
        *scratch_n = sources.n();

        char *scratch_ptr = scratch + sizeof(size_t) + sizeof(hpx_addr_t);
        for (size_t i = 0; i < sources.n(); ++i) {
          scratch_ptr = (char *)manager->serialize(&sref[i], scratch_ptr);
        }
      }

      int my_rank = hpx_get_my_rank();
      auto begin = parts->out_edges.begin();
      auto end = parts->out_edges.end();
      while (begin != end) {
        int curr_rank = begin->target->locality;
        auto curr = begin;
        while (curr != end && curr->target->locality == curr_rank) {
          ++curr;
        }

        //copy in edge data
        char *edgedata = scratch + header_size;
        size_t *edgecount = reinterpret_cast<size_t *>(edgedata);
        *edgecount = curr - begin;

        DAGInstigationRecord *edgerecords
            = reinterpret_cast<DAGInstigationRecord *>(edgedata
                                                        + sizeof(size_t));
        int i = 0;
        for (auto loop = begin; loop != curr; ++loop) {
          edgerecords[i].op = loop->op;
          edgerecords[i].target = loop->target->global_addx;
          edgerecords[i].idx = loop->target->index();
          ++i;
        }

        // Send parcel or do the work
        if (curr_rank == my_rank) {
          instigate_dag_eval_work(sources.n(), sref, tree->domain_,
                                  *edgecount, edgerecords);
        } else {
          size_t parcel_size = header_size + sizeof(size_t)
                               + sizeof(DAGInstigationRecord) * (*edgecount);
          hpx_parcel_t *parc = hpx_parcel_acquire(scratch, parcel_size);
          hpx_parcel_set_action(parc, instigate_dag_eval_remote_);
          hpx_parcel_set_target(parc, HPX_THERE(curr_rank));
          hpx_parcel_send_sync(parc);
        }

        begin = curr;
      }

      delete [] scratch;
    }

    return HPX_SUCCESS;
  }

  // TODO: Get this out of DualTree
  /// Action on remote side for DAG instigation
  ///
  /// Similar to the out edges for the Expansion LCOs, the edges out of the
  /// source nodes are sorted and the data is sent across the network once to
  /// each locality. This action handles the fan out once the data reaches
  /// the target locality.
  ///
  /// \param message - the message data
  /// \param bytes - the message size
  ///
  /// \returns - HPX_SUCCESS
  static int instigate_dag_eval_remote_handler(char *message, size_t bytes) {
    // unpack message into arguments to the local work function
    hpx_addr_t tree_addx = *(reinterpret_cast<hpx_addr_t *>(message));
    RankWise<dualtree_t> global_tree{tree_addx};
    auto local_tree = global_tree.here();

    size_t n_src = *(reinterpret_cast<size_t *>(message + sizeof(hpx_addr_t)));

    Source *sources = new Source[n_src];
    char *msg_ptr = message + sizeof(size_t) + sizeof(hpx_addr_t);
    {
      Array<Source> sarr{local_tree->source_gas};
      Serializer *manager = sarr.get_manager();
      for (size_t i = 0; i < n_src; ++i) {
        msg_ptr = (char *)manager->deserialize(msg_ptr, &sources[i]);
      }
    }

    size_t n_edges = *(reinterpret_cast<size_t *>(msg_ptr));
    DAGInstigationRecord *edges
        = reinterpret_cast<DAGInstigationRecord *>(msg_ptr + sizeof(size_t));

    // Detect if the edges have unknown target addresses and lookup the
    // correct address
    auto ttree = local_tree->target_tree_.here();
    for (size_t i = 0; i < n_edges; ++i) {
      if (edges[i].target == HPX_NULL) {
        edges[i].target = ttree->lookup_lco_addx(edges[i].idx, edges[i].op);
      }
    }

    instigate_dag_eval_work(n_src, sources, local_tree->domain_,
                            n_edges, edges);

    delete [] sources;

    return HPX_SUCCESS;
  }

  // TODO: Get this out of DualTree
  /// Perform the actual work of DAG instigation
  ///
  /// \param n_src - the number of sources
  /// \param sources - the source records
  /// \param domain - the domain geometry
  /// \param n_edges - the number of edges to process
  /// \param edge - the edge data
  static void instigate_dag_eval_work(size_t n_src,
                                      Source *sources,
                                      DomainGeometry &domain,
                                      size_t n_edges,
                                      DAGInstigationRecord *edge) {
    // loop over edges
    for (size_t i = 0; i < n_edges; ++i) {
      switch (edge[i].op) {
        case Operation::Nop:
          assert(0 && "Trouble handling DAG instigation");
          break;
        case Operation::StoM:
          {
            expansionlco_t expand{edge[i].target};
            Point center = domain.center_from_index(edge[i].idx);
            expand.S_to_M(center, sources, n_src, edge[i].idx);
          }
          break;
        case Operation::StoL:
          {
            expansionlco_t expand{edge[i].target};
            Point center = domain.center_from_index(edge[i].idx);
            expand.S_to_L(center, sources, n_src, edge[i].idx);
          }
          break;
        case Operation::MtoM:   // NOTE: Fall-through
        case Operation::MtoL:   //   |
        case Operation::LtoL:   //   |
        case Operation::MtoT:   //   |
        case Operation::LtoT:   //   v
          assert(0 && "Trouble handling DAG instigation");
          break;
        case Operation::StoT:
          {
            // S_to_T on expansion LCOs do not need any of the
            // expansionlco_t's state, so we create a default object.
            expansionlco_t expand{HPX_NULL};
            targetlco_t targets{edge[i].target};
            expand.S_to_T(sources, n_src, targets);
          }
          break;
        default:
          assert(0 && "Trouble handling DAG instigation");
          break;
      }
    }
  }

  // TODO: Get this out of DualTree
  /// Action to apply Method::aggregate
  ///
  /// This is an action that is called dependent on the children of this
  /// node all completing their Method application work.
  ///
  /// \param tree - the DualTree
  /// \param node - the node of the source tree
  /// \param done - the completion LCO that is set by this action
  ///
  /// \returns - HPX_SUCCESS
  static int source_apply_method_child_done_handler(dualtree_t *tree,
                                                    sourcenode_t *node,
                                                    hpx_addr_t done) {
    tree->method_.aggregate(node, &tree->domain_);
    int loc{0};
    if (node->idx.level() >= tree->unif_level_) {
      int dag_idx = sourcetree_t::get_unif_grid_index(node->idx,
                                                      tree->unif_level_);
      loc = tree->rank_of_unif_grid(dag_idx);
    }

    method_t::distropolicy_t::assign_for_source(node->dag, loc);

    hpx_lco_and_set(done, HPX_NULL);

    return HPX_SUCCESS;
  }

  // TODO: Get this out of DualTree
  /// Action to apply Method::generate
  ///
  /// This is a recursive spawn through the source tree.
  ///
  /// \param tree - the DualTree
  /// \param node - the node of the source tree
  /// \param done - LCO to signal when application of the method is complete
  ///
  /// \returns - HPX_SUCCESS
  static int source_apply_method_handler(dualtree_t *tree,
                                         sourcenode_t *node,
                                         hpx_addr_t done) {
    int n_children = node->n_children();
    if (n_children == 0) {
      tree->method_.generate(node, &tree->domain_);

      int dag_idx = sourcetree_t::get_unif_grid_index(node->idx,
                                                      tree->unif_level_);
      assert(dag_idx >= 0);
      int dag_rank = tree->rank_of_unif_grid(dag_idx);
      node->dag.set_parts_locality(dag_rank);
      node->dag.set_normal_locality(dag_rank);

      method_t::distropolicy_t::assign_for_source(node->dag, dag_rank);

      hpx_lco_and_set(done, HPX_NULL);

      return HPX_SUCCESS;
    }

    hpx_addr_t cdone = hpx_lco_and_new(n_children);
    assert(cdone != HPX_NULL);

    for (int i = 0; i < 8; ++i) {
      if (node->child[i] == nullptr) continue;
      hpx_call(HPX_HERE, source_apply_method_, HPX_NULL,
               &tree, &node->child[i], &cdone);
    }

    // Once the children are done, call aggregate here, deleting cdone as
    // the continuation.
    hpx_call_when_with_continuation(cdone,
                                    HPX_HERE, source_apply_method_child_done_,
                                    cdone, hpx_lco_delete_action,
                                    &tree, &node, &done);

    return HPX_SUCCESS;
  }

  // TODO: Get this out of DualTree
  /// Action to apply Method::inherit and Method::process
  ///
  /// This is a parallel spawn through the tree to apply the method to the
  /// given tree. This will generate the DAG.
  ///
  /// \param tree - the DualTree
  /// \param node - the target node in question
  /// \param consider - the list of source nodes to consider in process
  /// \param same_sandt - is S == T for this evaluation
  /// \param done - LCO to signal with completion
  ///
  /// \returns - HPX_SUCCESS
  static int target_apply_method_handler(dualtree_t *tree,
                                         targetnode_t *node,
                                         std::vector<sourcenode_t *> *consider,
                                         int same_sandt,
                                         hpx_addr_t done) {
    bool refine = false;
    if (node->idx.level() < tree->unif_level_) {
      refine = true;
    } else if (node->n_children()) {
      refine = tree->method_.refine_test((bool)same_sandt, node, *consider);
    }

    tree->method_.inherit(node, &tree->domain_, !refine);
    tree->method_.process(node, *consider, !refine, &tree->domain_);

    if (!refine) {
      // If we are not refining, then we are at a target leaf. This means
      // we should set the particle and normal DAG nodes to have this locality.

      int dag_idx = targettree_t::get_unif_grid_index(node->idx,
                                                      tree->unif_level_);
      int dag_rank = tree->rank_of_unif_grid(dag_idx);
      node->dag.set_parts_locality(dag_rank);
      node->dag.set_normal_locality(dag_rank);

      hpx_lco_set_lsync(done, 0, nullptr, HPX_NULL);
    } else {
      int n_children = node->n_children();
      hpx_addr_t cdone = hpx_lco_and_new(n_children);
      assert(cdone);

      for (int i = 0; i < 8; ++i) {
        if (node->child[i] == nullptr) continue;
        std::vector<sourcenode_t *> *ccons =
            new std::vector<sourcenode_t *>{};
        *ccons = *consider;
        hpx_call(HPX_HERE, target_apply_method_, HPX_NULL,
                 &tree, &node->child[i], &ccons, &same_sandt, &cdone);
      }

      int loc{0};
      if (node->idx.level() >= tree->unif_level_) {
        int dag_idx = targettree_t::get_unif_grid_index(node->idx,
                                                        tree->unif_level_);
        loc = tree->rank_of_unif_grid(dag_idx);
      }
      method_t::distropolicy_t::assign_for_target(node->dag, loc);

      assert(cdone != HPX_NULL);
      hpx_call_when(cdone, cdone, hpx_lco_delete_action,
                    done, nullptr, 0);
    }

    delete consider;

    return HPX_SUCCESS;
  }

  // TODO: Get this out of DualTree
  /// Action to set up termination detection for the DAG evaluation
  ///
  /// \param done - LCO for termination detection
  /// \param targs - target DAG nodes
  /// \param n_targs - number of target DAG nodes
  /// \param tint - internal target DAG nodes
  /// \param n_tint - number of internal target DAG nodes
  /// \param sint - internal source DAG nodes
  /// \param n_sint - number of internal source DAG nodes
  ///
  /// \returns - HPX_SUCCESS
  static int termination_detection_handler(hpx_addr_t done,
                                           std::vector<DAGNode *> *targs,
                                           size_t n_targs,
                                           std::vector<DAGNode *> *tint,
                                           size_t n_tint,
                                           std::vector<DAGNode *> *sint,
                                           size_t n_sint) {
    // If this is insufficiently parallel, we can always make this action
    // call itself with smaller and smaller chunks of the array.
    int nskip{0};

    int myrank = hpx_get_my_rank();
    for (size_t i = 0; i < n_targs; ++i) {
      assert((*targs)[i] != nullptr);
      if ((*targs)[i]->locality == myrank) {
        assert((*targs)[i]->global_addx != HPX_NULL);
        hpx_call_when((*targs)[i]->global_addx, done, hpx_lco_set_action,
                      HPX_NULL, nullptr, 0);
      } else {
        ++nskip;
      }
    }

    for (size_t i = 0; i < n_tint; ++i) {
      assert((*tint)[i] != nullptr);
      if ((*tint)[i]->locality == myrank) {
        assert((*tint)[i]->global_addx != HPX_NULL);
        hpx_call_when((*tint)[i]->global_addx, done, hpx_lco_set_action,
                      HPX_NULL, nullptr, 0);
      } else {
        ++nskip;
      }
    }

    for (size_t i = 0; i < n_sint; ++i) {
      assert((*sint)[i] != nullptr);
      if ((*sint)[i]->locality == myrank) {
        assert((*sint)[i]->global_addx != HPX_NULL);
        hpx_call_when((*sint)[i]->global_addx, done, hpx_lco_set_action,
                      HPX_NULL, nullptr, 0);
      } else {
        ++nskip;
      }
    }

    hpx_lco_and_set_num(done, nskip, HPX_NULL);

    return HPX_SUCCESS;
  }

  // TODO: Get this out of DualTree
  /// Action to destroy the DAG LCOs
  ///
  /// \param nodes - the DAG nodes
  /// \param n_nodes - the number of nodes
  /// \param type - nonzero for TargetLCO, zero for ExpansionLCO
  ///
  /// \returns - HPX_SUCCESS
  static int destroy_DAG_LCOs_handler(DAGNode **nodes,
                                      size_t n_nodes,
                                      int type) {
    int myrank = hpx_get_my_rank();
    for (size_t i = 0; i < n_nodes; ++i) {
      if (nodes[i]->locality == myrank) {
        assert(nodes[i]->global_addx != HPX_NULL);
        if (type) {
          // NOTE: destroy() does not care about the second argument to the
          // constructor, so we put in some junk
          auto temp = targetlco_t{nodes[i]->global_addx};
          temp.destroy();
        } else {
          auto temp = expansionlco_t{nodes[i]->global_addx};
          temp.destroy();
        }
      }
    }

    return HPX_SUCCESS;
  }

  //
  // Now for the data members
  //

  DomainGeometry domain_;     /// domain size
  int refinement_limit_;      /// refinement threshold
  int unif_level_;            /// level of uniform partition
  int dim3_;                  /// number of uniform nodes
  int same_sandt_;            /// Made from the same sources and targets
  hpx_addr_t unif_count_;     /// LCO reducing the uniform counts
  int *unif_count_value_;     /// local data storing the uniform counts
  int *rank_map_;             /// map unif grid index to rank
  method_t method_;           /// method used during DAG discovery

  RankWise<sourcetree_t> source_tree_;
  RankWise<targettree_t> target_tree_;

  hpx_addr_t source_gas;      /// the source array meta data
  hpx_addr_t target_gas;      /// the target array meta data


  static hpx_action_t domain_geometry_init_;
  static hpx_action_t domain_geometry_op_;
  static hpx_action_t set_domain_geometry_;
  static hpx_action_t init_partition_;
  static hpx_action_t recv_points_;
  static hpx_action_t send_points_;
  static hpx_action_t create_dual_tree_;
  static hpx_action_t finalize_partition_;
  static hpx_action_t source_apply_method_;
  static hpx_action_t source_apply_method_child_done_;
  static hpx_action_t target_apply_method_;
  static hpx_action_t destroy_DAG_LCOs_;
  static hpx_action_t termination_detection_;
  static hpx_action_t create_S_expansions_from_DAG_;
  static hpx_action_t create_T_expansions_from_DAG_;
  static hpx_action_t instigate_dag_eval_;
  static hpx_action_t instigate_dag_eval_remote_;
};

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::domain_geometry_init_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::domain_geometry_op_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::set_domain_geometry_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::init_partition_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::recv_points_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::send_points_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::create_dual_tree_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::finalize_partition_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::source_apply_method_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::source_apply_method_child_done_ =
    HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::target_apply_method_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::destroy_DAG_LCOs_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::termination_detection_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::create_S_expansions_from_DAG_ =
    HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::create_T_expansions_from_DAG_ =
    HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::instigate_dag_eval_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t DualTree<S, T, E, M>::instigate_dag_eval_remote_ =
    HPX_ACTION_NULL;


} // dashmm


#endif // __DASHMM_DUALTREE_H__
