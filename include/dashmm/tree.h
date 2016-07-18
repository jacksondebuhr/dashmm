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


#ifndef __DASHMM_TREE_H__
#define __DASHMM_TREE_H__


/// \file include/dashmm/tree.h
/// \brief Tree type


#include "dashmm/arrayref.h"
#include "dashmm/domaingeometry.h"
#include "dashmm/expansionlco.h"
#include "dashmm/dag.h"
#include "dashmm/shareddata.h"
#include "dashmm/targetlco.h"


namespace dashmm {


/// Forward declaration of Evaluator so that Tree can become friends
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class Method,
          typename DistroPolicy>
class Evaluator;


/// Node of either the source or target tree
///
/// This class is parameterized with the typical STEMD types for DASHMM, but
/// also includes one more parameter, which will always be either Source or
/// Target. The only different in the two trees is that the reference to the
/// array is either an array of Sources or Targets, hence the extra parameter.
///
/// There is very little in the way of behavior attached to this object, so it
/// exposes its data directly.
///
/// This object should only be created or destroyed from inside an HPX-5
/// thread.
template <typename Source, typename Target, typename Record,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class Method,
          typename DistroPolicy>
class TreeNode {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Method<Source, Target, Expansion, DistroPolicy>;

  using treenode_t = TreeNode<Source, Target, Record, Expansion, Method,
                              DistroPolicy>;
  using arrayref_t = ArrayRef<Record>;

  TreeNode(Index i, treenode_t *p)
      : idx{i}, parent{p}, child{}, dag{i}, parts{} { }

  ~TreeNode() {
    for (int i = 0; i < 8; ++i) {
      if (child[i] != nullptr) {
        delete child[i];
      }
    }
  }

  int n_children() const {
    int retval{0};
    for (int i = 0; i < 8; ++i) {
      if (child[i] != nullptr) {
        ++retval;
      }
    }
    return retval;
  }

  bool is_leaf() const {return n_children() == 0;}

  Index idx;                // Index of the node
  treenode_t *parent;       // The parent of this node
  treenode_t *child[8];     // The children of this node
  DAGInfo dag;              // The DAG info for this node
  arrayref_t parts;         // Reference to the particles (Source or Target)
};


/// Tree object for DASHMM evaluations
///
/// The tree object contains both the source and target trees, as well as some
/// information common to both structures. As this object creates the trees,
/// it will generate the lightweight DAG representation as well. This will be
/// useful for performing rapid work analyzing the DAG itself, without having
/// to fully instantiate the data that would inhabit the nodes of the DAG.
///
/// The methods in this object should only be called from HPX-5 threads. There
/// are some exceptions to this, but for the sake of simplicity, think of this
/// object as being useful only from HPX-5 threads. This does not preclude
/// the use of the address of such an object outside of HPX-5.
///
/// As the trees are created, the sources and targets will be sorted
/// geometrically consistent with the resulting tree.
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class Method,
          typename DistroPolicy>
class Tree {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Method<Source, Target, Expansion, DistroPolicy>;

  using sourcenode_t = TreeNode<Source, Target, Source, Expansion, Method,
                                DistroPolicy>;
  using targetnode_t = TreeNode<Source, Target, Target, Expansion, Method,
                                DistroPolicy>;
  using sourceref_t = ArrayRef<Source>;
  using targetref_t = ArrayRef<Target>;
  using tree_t = Tree<Source, Target, Expansion, Method, DistroPolicy>;
  using expansionlco_t = ExpansionLCO<Source, Target, Expansion, Method,
                                      DistroPolicy>;
  using targetlco_t = TargetLCO<Source, Target, Expansion, Method, DistroPolicy>;

  Tree(method_t met, size_t limit, int digits)
      : method_{met}, refinement_limit_{limit}, n_digits_{digits},
        source_root_{nullptr}, target_root_{nullptr}, domain_{nullptr} { }

  ~Tree() {
    if (source_root_ != nullptr) {
      delete source_root_;
    }
    if (target_root_ != nullptr) {
      delete target_root_;
    }
    domain_.destroy();
  }

  // delete the copy constructor and assignment
  Tree(const Tree &other) = delete;
  Tree &operator=(const Tree &other) = delete;

  // define move construction and assignment
  Tree(Tree &&other) {
    method_ = other.method;
    refinement_limit_ = other.refinement_limit_;
    n_digits_ = other.n_digits_;
    source_root_ = other.source_root_;
    target_root_ = other.target_root_;
    domain_ = other.domain_;

    other.source_root_ = nullptr;
    other.target_root_ = nullptr;
  }

  Tree &operator=(Tree &&other) {
    method_ = other.method;
    refinement_limit_ = other.refinement_limit_;
    n_digits_ = other.n_digits_;
    source_root_ = other.source_root_;
    target_root_ = other.target_root_;
    domain_ = other.domain_;

    other.source_root_ = nullptr;
    other.target_root_ = nullptr;

    return *this;
  }

  /// Returns the domain represented by the tree
  LocalData<DomainGeometry> domain() const {return domain_.value();}

  /// Returns the method to use for constructing the DAG
  const method_t &method() const {return method;}

  /// Returns the refinement limit for tree construction
  int refinement_limit() const {return refinement_limit_;}

  /// Returns the number of digits in any created expansion
  int n_digits() const {return n_digits_;}


  /// Compute the domain for the given sources and targets
  ///
  /// This will find the cube that surrounds the given source and target
  /// points. This sets the domain of this Tree. This routine must be
  /// called before the trees are created, otherwise undefined results will
  /// occur.
  ///
  /// \param sources - a reference to the source points
  /// \param targets - a reference to the target points
  /// \param same - indicates if the sources and targets are identical
  ///
  /// This has the side effect of setting this object's domain. The result
  /// can be retrieved using domain().
  void compute_domain(sourceref_t sources, targetref_t targets, bool same) {
    // Get source bounds
    hpx_addr_t srcbnd = hpx_lco_future_new(sizeof(BoundsResult));
    assert(srcbnd != HPX_NULL);
    hpx_addr_t sdata = sources.data();
    size_t scount = sources.n();
    hpx_call(HPX_THERE(0), source_bounds_, srcbnd, &sdata, &scount);

    BoundsResult bounds{Point{0.0, 0.0, 0.0}, Point{0.0, 0.0, 0.0}};
    hpx_lco_get(srcbnd, sizeof(BoundsResult), &bounds);
    hpx_lco_delete_sync(srcbnd);

    // get target bounds - if targets != sources
    if (!same) {
      hpx_addr_t trgbnd = hpx_lco_future_new(sizeof(BoundsResult));
      assert(trgbnd != HPX_NULL);
      hpx_addr_t tdata = targets.data();
      size_t tcount = targets.n();
      hpx_call(HPX_THERE(0), target_bounds_, trgbnd, &tdata, &tcount);

      BoundsResult otherbounds{Point{0.0, 0.0, 0.0}, Point{0.0, 0.0, 0.0}};
      hpx_lco_get(trgbnd, sizeof(BoundsResult), &otherbounds);
      bounds.low.lower_bound(otherbounds.low);
      bounds.high.upper_bound(otherbounds.high);

      hpx_lco_delete_sync(trgbnd);
    }

    // Finally, set the domain
    DomainGeometry dom = DomainGeometry{bounds.low, bounds.high, 1.0002};
    domain_.reset(&dom);
  }

  /// Partition the source tree
  ///
  /// Given the sources, this will create the hierarchical subdivision of the
  /// domain of this object. The domain of the tree should be set before this
  /// is called. Further, during the partitioning, the method will be applied
  /// and the portions of the DAG that connect in the source tree will be
  /// constructed.
  ///
  /// Please note that this will sort the source records geometrically to aid
  /// in the efficient evaluation of the multipole method.
  ///
  /// \param sources - the reference to the source particles
  ///
  /// \returns - an LCO that will trigger when the partitioning is complete.
  hpx_addr_t partition_source_tree(sourceref_t sources) {
    hpx_addr_t retval = hpx_lco_future_new(0);
    assert(retval != HPX_NULL);

    source_root_ = new sourcenode_t{Index{0, 0, 0, 0}, nullptr};

    SourcePartitionParams input{retval, this, source_root_, sources};

    hpx_call(HPX_THERE(0), source_partition_, HPX_NULL, &input, sizeof(input));

    return retval;
  }

  /// Partition the target tree
  ///
  /// The hierarchical subdivision of the targets is computed by this routine.
  /// Like with partition_source_tree(), this will sort the targets. The
  /// domain of the tree should have been computed before calling this function.
  /// As this partitioning is computed, the remainder of the DAG will be set up.
  ///
  /// \param targets - the reference to the target particles
  /// \param same_sources_and_targets - true if the sources and targets are
  ///                  identical.
  /// \param consider - a list of nodes to pass into the method::process
  ///                  function. Typically, this is just the source root.
  ///
  /// \returns - an LCO to indicate the completion of the partitioning
  // TODO remove consider from this
  hpx_addr_t partition_target_tree(targetref_t targets,
                                   bool same_sources_and_targets,
                                   std::vector<sourcenode_t *> *consider) {
    hpx_addr_t partdone = hpx_lco_future_new(0);
    assert(partdone != HPX_NULL);

    target_root_ = new targetnode_t{Index{0, 0, 0, 0}, nullptr};

    TargetPartitionParams input{};
    input.tree = this;
    input.node = target_root_;
    input.partdone = partdone;
    input.same_sources_and_targets = same_sources_and_targets;
    input.targets = targets;
    input.consider = consider;

    hpx_call(HPX_THERE(0), target_partition_, HPX_NULL, &input, sizeof(input));

    return partdone;
  }

  /// Partition the tree given the provided source and target points
  ///
  /// This will compute not only the bounding geometry of the points, but
  /// will also partition the tree. This operation is synchronous.
  ///
  /// \param sources - the reference to the source records
  /// \param targets - the reference to the target records
  void partition(sourceref_t sources, targetref_t targets, bool same_sandt) {
    compute_domain(sources, targets, same_sandt);

    hpx_addr_t sourcedone = partition_source_tree(sources);
    hpx_lco_wait(sourcedone);
    hpx_lco_delete_sync(sourcedone);

    std::vector<sourcenode_t *> *consider = new std::vector<sourcenode_t *>{};
    consider->push_back(source_root_);
    hpx_addr_t targetpartdone = partition_target_tree(targets, same_sandt,
                                                      consider);
    hpx_lco_wait(targetpartdone);
    hpx_lco_delete_sync(targetpartdone);
  }

  /// Create the DAG for this tree using the method specified at creation.
  ///
  /// This will allocate and collect the DAG nodes into the returned object.
  ///
  /// \returns - the resulting DAG.
  DAG create_DAG(bool same_sandt) {
    // Do work on the source tree
    hpx_addr_t sdone = hpx_lco_future_new(0);
    assert(sdone != HPX_NULL);

    tree_t *thetree = this;
    hpx_call(HPX_HERE, source_apply_method_, HPX_NULL,
             &thetree, &source_root_, &sdone);

    hpx_lco_wait(sdone);
    hpx_lco_delete_sync(sdone);

    // Do work on the target tree
    hpx_addr_t tdone = hpx_lco_future_new(0);
    assert(tdone != HPX_NULL);

    std::vector<sourcenode_t *> *consider = new std::vector<sourcenode_t *>{};
    consider->push_back(source_root_);
    int samearg = same_sandt ? 1 : 0;  // No HPX_BOOL type for action args
    hpx_call(HPX_HERE, target_apply_method_, HPX_NULL,
             &thetree, &target_root_, &consider, &samearg, &tdone);

    hpx_lco_wait(tdone);
    hpx_lco_delete_sync(tdone);

    DAG retval{};
    collect_DAG_nodes(retval);
    return retval;
  }

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
  void collect_DAG_nodes(DAG &dag) {
    collect_DAG_nodes_from_S_node(source_root_, dag.source_leaves,
                                  dag.source_nodes);
    collect_DAG_nodes_from_T_node(target_root_, dag.target_leaves,
                                  dag.target_nodes);
  }

  /// Create the LCOs from the DAG
  ///
  /// This will traverse the source and target tree creating any needed
  /// expansion LCOs. Further, it will create the target LCOs in the target
  /// tree.
  ///
  /// This is a synchronous operation.
  ///
  /// \param n_digits - the accuracy parameter of the expansion being used in
  ///                   this evaluation.
  void create_expansions_from_DAG(int n_digits) {
    hpx_addr_t done = hpx_lco_and_new(2);
    assert(done != HPX_NULL);

    tree_t *argthis = this;
    hpx_call(HPX_HERE, create_S_expansions_from_DAG_, HPX_NULL,
             &done, &n_digits, &argthis, &source_root_);
    hpx_call(HPX_HERE, create_T_expansions_from_DAG_, HPX_NULL,
             &done, &n_digits, &argthis, &target_root_);

    hpx_lco_wait(done);
    hpx_lco_delete_sync(done);
  }

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
  hpx_addr_t setup_termination_detection(DAG &dag) {
    size_t n_targs = dag.target_leaves.size();
    size_t n_tinternal = dag.target_nodes.size();
    size_t n_sinternal = dag.source_nodes.size();

    hpx_addr_t retval = hpx_lco_and_new(n_targs + n_tinternal + n_sinternal);
    assert(retval != HPX_NULL);

    std::vector<DAGNode *> *argaddx = &dag.target_leaves;
    std::vector<DAGNode *> *tinternalsaddx = &dag.target_nodes;
    std::vector<DAGNode *> *sinternalsaddx = &dag.source_nodes;
    hpx_call(HPX_HERE, termination_detection_, HPX_NULL, &retval,
             &argaddx, &n_targs, &tinternalsaddx, &n_tinternal,
             &sinternalsaddx, &n_sinternal);

    return retval;
  }

  /// Sets the edge lists of the expansion LCOs
  ///
  /// This is a separate phase as the DAG is constructed before the LCOs
  /// serving the work for a DAG node are created. Once the expansion LCOs
  /// are created, the addresses can then be shared with any LCO that needs
  /// that information.
  ///
  /// This is an asynchronous operation. The work will have started when this
  /// function returns, but it may not have ended. There is no returned LCO
  /// as the termination detection cannot trigger before this is done.
  ///
  /// \param dag - DAG object
  void setup_edge_lists(DAG &dag) {
    DAGNode **sdata = dag.source_nodes.data();
    size_t n_snodes = dag.source_nodes.size();
    DAGNode **tdata = dag.target_nodes.data();
    size_t n_tnodes = dag.target_nodes.size();
    hpx_call(HPX_HERE, edge_lists_, HPX_NULL,
             &sdata, &n_snodes, &tdata, &n_tnodes);
  }

  /// Initiate the DAG evaluation
  ///
  /// This starts the work of the evalution by starting the S->* work at the
  /// source nodes of the DAG.
  ///
  /// This is an asynchronous operation. The termination detection cannot
  /// possibly trigger before this is completed, so waiting on the termination
  /// of the full evaluation implicitly waits on this operation.
  void start_DAG_evaluation() {
    tree_t *targ = this;
    hpx_call(HPX_HERE, instigate_dag_eval_, HPX_NULL, &targ, &source_root_);
  }

  /// Destroys the LCOs associated with the DAG
  ///
  /// This is a synchronous operation. This destroys not only the expansion
  /// LCOs, but also the target LCOs.
  ///
  /// \param targets - the target nodes of the DAG
  /// \param internal - the internal nodes of the DAG
  void destroy_DAG_LCOs(DAG &dag) {
    hpx_addr_t done = hpx_lco_and_new(3);
    assert(done != HPX_NULL);

    DAGNode **data = dag.target_leaves.data();
    size_t n_data = dag.target_leaves.size();
    hpx_call(HPX_HERE, destroy_DAG_LCOs_, done, &data, &n_data);

    data = dag.target_nodes.data();
    n_data = dag.target_nodes.size();
    hpx_call(HPX_HERE, destroy_DAG_LCOs_, done, &data, &n_data);

    data = dag.source_nodes.data();
    n_data = dag.source_nodes.size();
    hpx_call(HPX_HERE, destroy_DAG_LCOs_, done, &data, &n_data);

    hpx_lco_wait(done);
    hpx_lco_delete_sync(done);
  }

 private:
  friend class Evaluator<Source, Target, Expansion, Method, DistroPolicy>;

  /// The result of finding bounds
  struct BoundsResult {
    Point low;
    Point high;
  };

  /// Marshalled argument type for source tree partitioning
  struct SourcePartitionParams {
    hpx_addr_t partdone;
    tree_t *tree;
    sourcenode_t *node;
    sourceref_t sources;
  };

  /// Marshalled argument type for target tree partitioning
  struct TargetPartitionParams {
    tree_t *tree;
    targetnode_t *node;
    hpx_addr_t partdone;
    bool same_sources_and_targets;
    targetref_t targets;
    // We are always working on one locality for the tree construction, so we
    // can get away with passing these pointers around. If we ever build the
    // tree across localities, this will need to change.
    std::vector<sourcenode_t *> *consider;
  };

  static int source_bounds_handler(hpx_addr_t data, size_t count) {
    source_t *user{nullptr};
    assert(hpx_gas_try_pin(data, (void **)&user));

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
    assert(hpx_gas_try_pin(data, (void **)&user));

    BoundsResult retval{Point{1.0e34, 1.0e34, 1.0e34},
                        Point{-1.0e34, -1.0e34, -1.0e34}};

    for (size_t i = 0; i < count; ++i) {
      retval.low.lower_bound(user[i].position);
      retval.high.upper_bound(user[i].position);
    }

    return HPX_THREAD_CONTINUE(retval);
  }

  static int source_partition_handler(SourcePartitionParams *parms,
                                      size_t UNUSED) {
    tree_t *tree = parms->tree;
    sourcenode_t *node = parms->node;

    auto domain = tree->domain_.value();

    assert(node->parts.data() == HPX_NULL);
    node->parts = parms->sources;

    // If we are at a leaf, we can generate
    if (parms->sources.n() <= tree->refinement_limit_) {
      hpx_lco_set(parms->partdone, 0, nullptr, HPX_NULL, HPX_NULL);
      return HPX_SUCCESS;
    }

    // Otherwise, partition the sources and spawn more work
    source_t *source_parts{nullptr};
    assert(hpx_gas_try_pin(parms->sources.data(), (void **)&source_parts));

    source_t *splits[9] { };
    splits[0] = source_parts;
    splits[8] = &source_parts[parms->sources.n()];

    LocalData<DomainGeometry> ldom = tree->domain_.value();
    Point cen{ldom->center_from_index(node->idx)};
    double z_center = cen.z();
    auto z_comp = [&z_center](source_t &a) {
                    return a.position.z() < z_center;
                  };
    splits[4] = std::partition(splits[0], splits[8], z_comp);

    double y_center = cen.y();
    auto y_comp = [&y_center](source_t &a) {
                    return a.position.y() < y_center;
                  };
    splits[2] = std::partition(splits[0], splits[4], y_comp);
    splits[6] = std::partition(splits[4], splits[8], y_comp);

    double x_center = cen.x();
    auto x_comp = [&x_center](source_t &a) {
                    return a.position.x() < x_center;
                  };
    splits[1] = std::partition(splits[0], splits[2], x_comp);
    splits[3] = std::partition(splits[2], splits[4], x_comp);
    splits[5] = std::partition(splits[4], splits[6], x_comp);
    splits[7] = std::partition(splits[6], splits[8], x_comp);

    hpx_gas_unpin(parms->sources.data());

    // Find some counts
    sourceref_t cparts[8] { };
    int n_children{0};
    {
      int n_offset{0};
      for (int i = 0; i < 8; ++i) {
        int n_per_child = splits[i + 1] - splits[i];
        if (n_per_child) {
          ++n_children;
          cparts[i] = parms->sources.slice(n_offset, n_per_child);
        }
        n_offset += n_per_child;
      }
    }

    hpx_addr_t childpartdone = hpx_lco_and_new(n_children);
    assert(childpartdone != HPX_NULL);

    // Many arguments to the children's partition are the same
    SourcePartitionParams args{ };
    args.partdone = childpartdone;
    args.tree = tree;
    for (int i = 0; i < 8; ++i) {
      if (!cparts[i].valid()) {
        node->child[i] = nullptr;
        continue;
      }

      node->child[i] = new sourcenode_t{node->idx.child(i), node};

      // We only set the stuff that changes here
      args.node = node->child[i];
      args.sources = cparts[i];
      hpx_call(HPX_HERE, source_partition_, HPX_NULL, &args, sizeof(args));
    }

    hpx_call_when(childpartdone, childpartdone, hpx_lco_delete_action,
                  parms->partdone, nullptr, 0);

    return HPX_SUCCESS;
  }

  static int target_partition_handler(TargetPartitionParams *parms,
                                      size_t UNUSED) {
    targetnode_t *node = parms->node;
    tree_t *tree = parms->tree;

    node->parts = parms->targets;

    if (parms->targets.n() > tree->refinement_limit_) {
      // partition
      target_t *T{nullptr};
      assert(hpx_gas_try_pin(parms->targets.data(), (void **)&T));
      target_t *splits[9] { };
      splits[0] = T;
      splits[8] = &T[parms->targets.n()];

      LocalData<DomainGeometry> ldom = tree->domain_.value();
      Point cen{ldom->center_from_index(node->idx)};
      double z_center = cen.z();
      auto z_comp = [&z_center](target_t &a) {
                      return a.position.z() < z_center;
                    };
      double y_center = cen.y();
      auto y_comp = [&y_center](target_t &a) {
                      return a.position.y() < y_center;
                    };
      double x_center = cen.x();
      auto x_comp = [&x_center](target_t &a) {
                      return a.position.x() < x_center;
                    };

      if (parms->same_sources_and_targets) {
        // NOTE: The passed in value of same_sources_and_targets had better
        // be correct, or there will be trouble.
        splits[4] = std::partition_point(splits[0], splits[8], z_comp);

        splits[2] = std::partition_point(splits[0], splits[4], y_comp);
        splits[6] = std::partition_point(splits[4], splits[8], y_comp);

        splits[1] = std::partition_point(splits[0], splits[2], x_comp);
        splits[3] = std::partition_point(splits[2], splits[4], x_comp);
        splits[5] = std::partition_point(splits[4], splits[6], x_comp);
        splits[7] = std::partition_point(splits[6], splits[8], x_comp);
      } else {
        splits[4] = std::partition(splits[0], splits[8], z_comp);

        splits[2] = std::partition(splits[0], splits[4], y_comp);
        splits[6] = std::partition(splits[4], splits[8], y_comp);

        splits[1] = std::partition(splits[0], splits[2], x_comp);
        splits[3] = std::partition(splits[2], splits[4], x_comp);
        splits[5] = std::partition(splits[4], splits[6], x_comp);
        splits[7] = std::partition(splits[6], splits[8], x_comp);
      }

      hpx_gas_unpin(parms->targets.data());

      targetref_t cparts[8] { };
      int n_children{0};
      {
        int n_offset{0};
        for (int i = 0; i < 8; ++i) {
          int n_per_child = splits[i + 1] - splits[i];
          if (n_per_child) {
            ++n_children;
            cparts[i] = parms->targets.slice(n_offset, n_per_child);
          }
          n_offset += n_per_child;
        }
      }

      hpx_addr_t partdone = hpx_lco_and_new(n_children);
      assert(partdone != HPX_NULL);

      // set up the arguments to the partition actions -- the constant parts
      TargetPartitionParams args{};
      args.tree = tree;
      args.partdone = partdone;
      args.same_sources_and_targets = parms->same_sources_and_targets;
      size_t argssize = sizeof(args);

      for (int i = 0; i < 8; ++i) {
        if (!cparts[i].valid()) {
          node->child[i] = nullptr;
          continue;
        }

        node->child[i] = new targetnode_t{node->idx.child(i), node};

        args.node = node->child[i];
        args.targets = cparts[i];
        args.consider = new std::vector<sourcenode_t *>{};
        *args.consider = *parms->consider; // TODO consider not needed

        hpx_call(HPX_HERE, target_partition_, HPX_NULL, &args, argssize);
      }

      // This propagates the fact that these partition actions are done back
      // up the tree.
      hpx_call_when(partdone, partdone, hpx_lco_delete_action,
                    parms->partdone, nullptr, 0);
    } else {
      hpx_lco_set_lsync(parms->partdone, 0, nullptr, HPX_NULL);
    }

    // Now we clear out the consider vector's memory
    delete parms->consider;

    return HPX_SUCCESS;
  }

  static int source_apply_method_child_done_handler(tree_t *tree,
                                                    sourcenode_t *node,
                                                    hpx_addr_t done) {
    auto domain = tree->domain_.value();
    tree->method_.aggregate(node, domain.value());
    hpx_lco_delete_sync(done);
    return HPX_SUCCESS;
  }

  static int source_apply_method_handler(tree_t *tree, sourcenode_t *node,
                                         hpx_addr_t done) {
    int n_children = node->n_children();
    if (n_children == 0) {
      auto domain = tree->domain_.value();
      tree->method_.generate(node, domain.value());
      node->dag.set_parts_locality(hpx_get_my_rank());
      hpx_lco_set(done, 0, nullptr, HPX_NULL, HPX_NULL);
      return HPX_SUCCESS;
    }

    hpx_addr_t cdone = hpx_lco_and_new(n_children);
    assert(cdone != HPX_NULL);

    for (int i = 0; i < 8; ++i) {
      if (node->child[i] == nullptr) continue;
      hpx_call(HPX_HERE, source_apply_method_, HPX_NULL,
               &tree, &node->child[i], &cdone);
    }

    // Once the children are done, call aggregate here, continuing a set to
    // done once that has happened.
    hpx_call_when(cdone, HPX_HERE, source_apply_method_child_done_, done,
                  &tree, &node, &cdone);

    return HPX_SUCCESS;
  }

  static int target_apply_method_handler(tree_t *tree, targetnode_t *node,
                                         std::vector<sourcenode_t *> *consider,
                                         int same_sandt, hpx_addr_t done) {
    bool refine = false;
    if (node->parts.n() > tree->refinement_limit_) {
      refine = tree->method_.refine_test((bool)same_sandt, node, *consider);
    }

    auto domain = tree->domain_.value();
    tree->method_.inherit(node, domain.value(), !refine);
    tree->method_.process(node, *consider, !refine, domain.value());
    node->dag.set_parts_locality(hpx_get_my_rank());

    if (refine) {
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

      hpx_call_when(cdone, cdone, hpx_lco_delete_action,
                    done, nullptr, 0);
    } else {
      hpx_lco_set_lsync(done, 0, nullptr, HPX_NULL);
    }

    delete consider;

    return HPX_SUCCESS;
  }

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

  static int create_S_expansions_from_DAG_handler(
        hpx_addr_t done, int n_digits, tree_t *tree, sourcenode_t *node) {
    auto domain = tree->domain_.value();
    Point n_center = domain->center_from_index(node->idx);

    // create the normal expansion if needed
    if (node->dag.has_normal()) {
      std::unique_ptr<expansion_t> input_expand{
        new expansion_t{n_center, n_digits, kSourcePrimary}
      };
      expansionlco_t expand(node->dag.normal()->in_edges.size(),
                            node->dag.normal()->out_edges.size(),
                            tree->domain_, node->idx, std::move(input_expand),
                            HPX_THERE(node->dag.normal()->locality));
      node->dag.set_normal_expansion(expand.lco(), expand.accuracy());
    }

    // If there is to be an intermediate expansion, create that
    if (node->dag.has_interm()) {
      std::unique_ptr<expansion_t> interm_expand{
        new expansion_t{n_center, n_digits, kSourceIntermediate}
      };
      expansionlco_t intexp_lco(node->dag.interm()->in_edges.size(),
                                node->dag.interm()->out_edges.size(),
                                tree->domain_, node->idx,
                                std::move(interm_expand),
                                HPX_THERE(node->dag.interm()->locality));
      node->dag.set_interm_expansion(intexp_lco.lco(), intexp_lco.accuracy());
    }

    // spawn work at children
    int n_children{node->n_children()};

    if (n_children) {
      hpx_addr_t cdone = hpx_lco_and_new(n_children);
      assert(cdone != HPX_NULL);

      for (int i = 0; i < 8; ++i) {
        if (node->child[i] != nullptr) {
          hpx_call(HPX_HERE, create_S_expansions_from_DAG_, HPX_NULL,
                   &cdone, &n_digits, &tree, &node->child[i]);
        }
      }

      // This will set the parent's LCO as well as delete cdone
      hpx_call_when(cdone, cdone, hpx_lco_delete_action,
                    done, nullptr, 0);
    } else {
      node->dag.set_sourceref(node->parts.data(), node->parts.n());

      hpx_lco_set(done, 0, nullptr, HPX_NULL, HPX_NULL);
    }

    return HPX_SUCCESS;
  }

  static int create_T_expansions_from_DAG_handler(
        hpx_addr_t done, int n_digits, tree_t *tree, targetnode_t *node) {
    auto domain = tree->domain_.value();
    Point n_center = domain->center_from_index(node->idx);

    // create the normal expansion if needed
    if (node->dag.has_normal()) {
      std::unique_ptr<expansion_t> input_expand{
        new expansion_t{n_center, n_digits, kTargetPrimary}
      };
      expansionlco_t expand(node->dag.normal()->in_edges.size(),
                            node->dag.normal()->out_edges.size(),
                            tree->domain_, node->idx, std::move(input_expand),
                            HPX_THERE(node->dag.normal()->locality));
      node->dag.set_normal_expansion(expand.lco(), expand.accuracy());
    }

    // If there is to be an intermediate expansion, create that
    if (node->dag.has_interm()) {
      std::unique_ptr<expansion_t> interm_expand{
        new expansion_t{n_center, n_digits, kTargetIntermediate}
      };
      expansionlco_t intexp_lco(node->dag.interm()->in_edges.size(),
                                node->dag.interm()->out_edges.size(),
                                tree->domain_, node->idx,
                                std::move(interm_expand),
                                HPX_THERE(node->dag.interm()->locality));
      node->dag.set_interm_expansion(intexp_lco.lco(), intexp_lco.accuracy());
    }

    // NOTE: this spawn through the tree does not end when the tree ends.
    // Instead, we have to check if this node has a parts node in the DAG.
    // If so, this branch is done, and we need not spawn more.

    // Here is where we make the target lco if needed
    if (node->dag.has_parts()) {
      targetlco_t tlco{node->dag.parts()->in_edges.size(), node->parts,
                       HPX_THERE(node->dag.parts()->locality)};
      node->dag.set_targetlco(tlco.lco(), tlco.n());

      hpx_lco_set(done, 0, nullptr, HPX_NULL, HPX_NULL);
    } else {
      hpx_addr_t cdone = hpx_lco_and_new(node->n_children());
      assert(cdone != HPX_NULL);

      for (int i = 0; i < 8; ++i) {
        if (node->child[i] != nullptr) {
          hpx_call(HPX_HERE, create_T_expansions_from_DAG_, HPX_NULL,
                   &cdone, &n_digits, &tree, &node->child[i]);
        }
      }

      // This will set the parent's LCO as well as delete cdone
      hpx_call_when(cdone, cdone, hpx_lco_delete_action,
                    done, nullptr, 0);
    }

    return HPX_SUCCESS;
  }

  static int edge_lists_handler(DAGNode **snodes, size_t n_snodes,
                                DAGNode **tnodes, size_t n_tnodes) {
    // If this is a bottleneck, we can easily make this some kind of parfor
    for (size_t i = 0; i < n_snodes; ++i) {
      expansionlco_t expand{snodes[i]->global_addx,
                            (int)snodes[i]->other_member};
      expand.set_out_edge_data(snodes[i]->out_edges);
    }
    for (size_t i = 0; i < n_tnodes; ++i) {
      expansionlco_t expand{tnodes[i]->global_addx,
                            (int)tnodes[i]->other_member};
      expand.set_out_edge_data(tnodes[i]->out_edges);
    }
    return HPX_SUCCESS;
  }

  static int instigate_dag_eval_handler(tree_t *tree, sourcenode_t *node) {
    int n_children{node->n_children()};
    if (n_children > 0) {
      for (int i = 0; i < 8; ++i) {
        if (node->child[i] != nullptr) {
          hpx_call(HPX_HERE, instigate_dag_eval_, HPX_NULL,
                  &tree, &node->child[i]);
        }
      }
    } else {
      // At a leaf, we do actual work
      DAGNode *parts = node->dag.parts();
      assert(parts != nullptr);
      sourceref_t sources = node->parts;

      // We first sort the out edges by locality
      std::sort(parts->out_edges.begin(), parts->out_edges.end(),
                DAG::compare_edge_locality);

      // loop over edges
      for (size_t i = 0; i < parts->out_edges.size(); ++i) {
        switch (parts->out_edges[i].op) {
          case Operation::Nop:
            assert(0 && "Trouble handling DAG instigation");
            break;
          case Operation::StoM:
            {
              const DAGNode *normal = node->dag.normal();
              assert(normal);
              expansionlco_t expand{normal->global_addx,
                                    (int)normal->other_member};
              auto domain = tree->domain_.value();
              double scale = 1.0 / domain->size_from_level(normal->idx.level());
              Point center = domain->center_from_index(normal->idx);
              expand.S_to_M(center, sources, scale);
            }
            break;
          case Operation::StoL:
            {
              // Get target DAGNode
              const DAGNode *normal = parts->out_edges[i].target;
              assert(normal);
              expansionlco_t expand{normal->global_addx,
                                    (int)normal->other_member};
              auto domain = tree->domain_.value();
              double scale = 1.0 / domain->size_from_level(normal->idx.level());
              Point center = domain->center_from_index(normal->idx);
              expand.S_to_L(center, sources, scale);
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
              expansionlco_t expand{HPX_NULL, 0};
              targetlco_t targets{parts->out_edges[i].target->global_addx,
                                  parts->out_edges[i].target->other_member};
              expand.S_to_T(sources, targets);
            }
            break;
          default:
            assert(0 && "Trouble handling DAG instigation");
            break;
        }
      }
    }

    return HPX_SUCCESS;
  }

  static int termination_detection_handler(hpx_addr_t done,
                                           std::vector<DAGNode *> *targs,
                                           size_t n_targs,
                                           std::vector<DAGNode *> *tint,
                                           size_t n_tint,
                                           std::vector<DAGNode *> *sint,
                                           size_t n_sint) {
    // If this is insufficiently parallel, we can always make this action
    // call itself with smaller and smaller chunks of the array.
    for (size_t i = 0; i < n_targs; ++i) {
      assert((*targs)[i] != nullptr);
      assert((*targs)[i]->global_addx != HPX_NULL);
      hpx_call_when((*targs)[i]->global_addx, done, hpx_lco_set_action,
                    HPX_NULL, nullptr, 0);
    }

    for (size_t i = 0; i < n_tint; ++i) {
      assert((*tint)[i] != nullptr);
      assert((*tint)[i]->global_addx != HPX_NULL);
      hpx_call_when((*tint)[i]->global_addx, done, hpx_lco_set_action,
                    HPX_NULL, nullptr, 0);
    }

    for (size_t i = 0; i < n_sint; ++i) {
      assert((*sint)[i] != nullptr);
      assert((*sint)[i]->global_addx != HPX_NULL);
      hpx_call_when((*sint)[i]->global_addx, done, hpx_lco_set_action,
                    HPX_NULL, nullptr, 0);
    }

    return HPX_SUCCESS;
  }

  static int destroy_DAG_LCOs_handler(DAGNode **nodes, size_t n_nodes) {
    // We could add more parallelism here if needed.

    // Some of these might be remote, so we use async delete on these LCOs.
    hpx_addr_t alldel = hpx_lco_and_new(n_nodes);
    assert(alldel != HPX_NULL);

    for (size_t i = 0; i < n_nodes; ++i) {
      assert(nodes[i]->global_addx != HPX_NULL);
      hpx_lco_delete(nodes[i]->global_addx, alldel);
    }

    hpx_lco_wait(alldel);
    hpx_lco_delete_sync(alldel);

    return HPX_SUCCESS;
  }

  // Data that is constant for each node of the tree
  method_t method_;
  size_t refinement_limit_;
  int n_digits_;

  // The roots of the two trees
  sourcenode_t *source_root_;
  targetnode_t *target_root_;

  // Members in SharedData
  SharedData<DomainGeometry> domain_;

  // Actions
  static hpx_action_t source_bounds_;
  static hpx_action_t target_bounds_;
  static hpx_action_t source_partition_;
  static hpx_action_t target_partition_;

  static hpx_action_t source_apply_method_child_done_;
  static hpx_action_t source_apply_method_;
  static hpx_action_t target_apply_method_;

  static hpx_action_t create_S_expansions_from_DAG_;
  static hpx_action_t create_T_expansions_from_DAG_;
  static hpx_action_t edge_lists_;
  static hpx_action_t instigate_dag_eval_;
  static hpx_action_t termination_detection_;
  static hpx_action_t destroy_DAG_LCOs_;
};


template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, E, M, D>::source_bounds_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, E, M, D>::target_bounds_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, E, M, D>::source_partition_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, E, M, D>::target_partition_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, E, M, D>::source_apply_method_child_done_ =
    HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, E, M, D>::source_apply_method_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, E, M, D>::target_apply_method_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, E, M, D>::create_S_expansions_from_DAG_ =
    HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, E, M, D>::create_T_expansions_from_DAG_ =
    HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, E, M, D>::edge_lists_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, E, M, D>::instigate_dag_eval_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, E, M, D>::termination_detection_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class,
                    typename> class M,
          typename D>
hpx_action_t Tree<S, T, E, M, D>::destroy_DAG_LCOs_ = HPX_ACTION_NULL;


} // namespace dashmm


#endif // __DASHMM_TREE_H__
