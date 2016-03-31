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


#include "dashmm/domaingeometry.h"
#include "dashmm/sourcenode.h"
#include "dashmm/targetnode.h"
#include "dashmm/treeshared.h"


namespace dashmm {


/// Forward declaration of Evaluator so that Tree can become friends
template <typename Source, typename Target,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class> class Method>
class Evaluator;


/// Node of either the source or target tree
///
/// This class is parameterized with the typical STEM types for DASHMM, but
/// also includes one more parameter, which will always be either Source or
/// Target. The only different in the two trees is that the reference to the
/// array is either an array of Sources or Targets, hence the extra parameter.
///
/// There is very little in the way of behavior attached to this object, so it
/// exposes its data directly.
///
/// NOTE: This object should only be created or destroyed from inside an HPX-5
/// thread.
template <typename Source, typename Target, typename Record,
          template <typename, typename> class Expansion,
          template <typename, typename,
                    template <typename, typename> class> class Method>
class TreeNode {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Method<Source, Target, Expansion>;

  using treenode_t = TreeNode<Source, Target, Record, Expansion, Method>;
  using arrayref_t = ArrayRef<Record>;

  TreeNode(Index i, treenode_t *p)
      : idx{i}, parent{p}, child{}, dag{}, parts{} { }

  ~TreeNode() {
    // TODO Do we want to do the following in parallel? That is, spawn work
    // for the children?
    for (int i = 0; i < 8; ++i) {
      if (child[i] != nullptr) {
        delete child[i];
      }
    }
  }

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
                    template <typename, typename> class> class Method>
class Tree {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = Expansion<Source, Target>;
  using method_t = Method<Source, Target, Expansion>;

  using sourcenode_t = TreeNode<Source, Target, Source, Expansion, Method>;
  using targetnode_t = TreeNode<Source, Target, Target, Expansion, Method>;
  using sourceref_t = ArrayRef<Source>;
  using targetref_t = ArrayRef<Target>;

  Tree(DomainGeometry dom, method_t met, int limit, int digits)
      : domain_{dom}, method_{met}, refinement_limit_{limit},
        n_digits_{digits}, source_root_{nullptr}, target_root_{nullptr} { }

  ~Tree() {
    if (source_root_ != nullptr) {
      delete source_root_;
    }
    if (target_root_ != nullptr) {
      delete target_root_;
    }
  }

  /// Returns the domain represented by the tree
  const DomainGeometry &domain() const {return domain_;}

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
    if (same) {
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
    domain_ = DomainGeometry{bounds.low, bounds.high, 1.0002};
  }

  // Document this - this also sets up the DAG stuff, or starts to
  hpx_addr_t partition_source_tree(sourceref_t sources) {
    hpx_addr_t retval = hpx_lco_future_new(0);
    assert(retval != HPX_NULL);

    source_root_ = new sourceroot_t{Index{0, 0, 0, 0}, nullptr};

    SourcePartitionParams input{retval, this, source_root_, sources};

    hpx_call(HPX_THERE(0), source_partition_, HPX_NULL, &input, sizeof(input));

    return retval;
  }

  // NOTE: Same as previous
  hpx_addr_t partition_target_tree(targetref_t targets,
                                   bool same_sources_and_targets,
                                   std::vector<sourcenode_t *> consider) {
    hpx_addr_t partdone = hpx_lco_future_new(0);
    assert(partdone != HPX_NULL);

    target_root_ = new targetnode_t{Index{0, 0, 0, 0}, nullptr};

    TargetPartitionParams input{};
    input.tree = this;
    input.node = target_root_;
    input.partdone = partdone;
    input.same_sources_and_targets = same_sources_and_targets;
    input.targets = targets;
    input.which_child = 0;
    input.consider = consider;

    hpx_call(HPX_THERE(0), target_partition_, HPX_NULL, &input, sizeof(input));

    return partdone;
  }

  // TODO: document
  void partition(sourceref_t sources, targetref_t targets) {
    bool same_sandt = (sources.data() == targets.data()
                        && sources.n() == targets.n());

    compute_domain(sources, targets, same_sandt);

    hpx_addr_t sourcedone = partition_source_tree(source);
    hpx_lco_wait(sourcedone);
    hpx_lco_delete_sync(sourcedone);

    std::vector<sourcenode_t *> consider{source_root_};
    hpx_addr_t targetpartdone = partition_target_tree(targets, same_sandt,
                                                      consider);
    hpx_lco_wait(targetpartdone);
    hpx_lco_delete_sync(targetpartdone);
  }

 private:
  friend class Evaluator<Source, Target, Expansion, Method>;

  /// The result of finding bounds
  struct BoundsResult {
    Point low;
    Point high;
  };

  struct SourcePartitionParams {
    hpx_addr_t partdone;
    tree_t *tree;
    sourcenode_t *node;
    sourceref_t sources;
  };

  struct TargetPartitionParams {
    tree_t *tree;
    targetnode_t *node;
    hpx_addr_t partdone;
    bool same_sources_and_targets;
    targetref_t targets;
    int which_child;
    // NOTE: Normally this is not sound at all. But the memory associated
    // with the actual vector elements will be local, as we are always
    // partitioning the tree on one locality. So this is why we can get away
    // with this for this action. If this is ever not done locally, then we
    // shall need to change this.
    std::vector<sourcenode_t *> consider;
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

  static int source_child_done_handler(tree_t *tree, sourcenode_t *node,
                                       hpx_addr_t partdone) {
    tree->method_.aggregate(node, tree->n_digits_);
    hpx_lco_delete_sync(partdone);
    return HPX_SUCCESS;
  }

  static int source_partition_handler(SourcePartitionParams *parms,
                                      size_t UNUSED) {
    tree_t *tree = parms->tree;
    sourcenode_t *node = parms->node;

    // If we are at a leaf, we can generate
    if (parms->sources.n() <= tree->refinement_limit_) {
      assert(node->sources.data() == HPX_NULL);
      node->sources = parms->sources;

      node->dag.add_parts();
      tree->method_.generate(node, tree->n_digits());

      hpx_lco_set(parms->partdone, 0, nullptr, HPX_NULL, HPX_NULL);

      return HPX_SUCCESS;
    }

    // Otherwise, partition the sources and spawn more work
    source_t *source_parts{nullptr};
    assert(hpx_gas_try_pin(parms->sources.data(), (void **)&source_parts));

    source_t *splits[9] { };
    splits[0] = source_parts;
    splits[8] = &source_parts[parms->sources.n()];

    Point cen{tree->domain_.center_from_index(node->idx)};
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
    PartitionParams args{ };
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

    // When the child has finished partitioning, we aggregate on this node
    hpx_call_when(childpartdone, HPX_HERE, source_child_done_,
                  parms->partdone, &tree, &node, &childpartdone);

    return HPX_SUCCESS;
  }

  static int target_partition_handler(TargetPartitionParams *parms,
                                      size_t UNUSED) {
    targetnode_t *node = parms->node;
    tree_t *tree = parms->tree;

    bool refine = false;
    if (parms->targets.n() > parms->limit) {
      refine = tree->method_.refine_test(parms->same_sources_and_targets, node,
                                         parms->consider);
    }

    if (!refine) {
      node->dag.add_parts();
    }

    tree->method_.inherit(node, tree->n_digits_, parms->which_child);
    tree->method_.process(node, parms->consider, !refine);

    if (refine) {
      // partition
      target_t *T{nullptr};
      assert(hpx_gas_try_pin(parms->targets.data(), (void **)&T));
      target_t *splits[9] { };
      splits[0] = T;
      splits[8] = &T[parms->targets.n()];

      Point cen{tree->domain_.center_from_index(node->idx)};
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
      args.consider = parms->consider;

      for (int i = 0; i < 8; ++i) {
        if (!cparts[i].valid()) {
          node->child[i] = nullptr;
          continue;
        }

        node->child[i] = new targetnode_t{node->idx.child(i), node};

        args.node = node->child[i];
        args.targets = cparts[i];
        args.which_child = i;

        hpx_call(HPX_HERE, target_partition_, HPX_NULL, &args, argssize);
      }
      delete [] args;

      // This propagates the fact that these partition actions are done back
      // up the tree.
      hpx_call_when(partdone, partdone, hpx_lco_delete_action,
                    parms->partdone, nullptr, 0);
    }

    if (!refine) {
      hpx_lco_set_lsync(parms->partdone, 0, nullptr, HPX_NULL);
    }

    return HPX_SUCCESS;
  }


  // Data that is constant for each node of the tree
  DomainGeometry domain_;
  method_t method_;
  int refinement_limit_;
  int n_digits_;

  // The roots of the two trees
  sourcenode_t *source_root_;
  targetnode_t *target_root_;

  // Actions
  static hpx_action_t source_bounds_;
  static hpx_action_t target_bounds_;
  static hpx_action_t source_child_done_;
  static hpx_action_t source_partition_;
  static hpx_action_t target_partition_;
};


template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Tree<S, T, E, M>::source_bounds_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Tree<S, T, E, M>::target_bounds_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Tree<S, T, E, M>::source_child_done_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Tree<S, T, E, M>::source_partition_ = HPX_ACTION_NULL;

template <typename S, typename T,
          template <typename, typename> class E,
          template <typename, typename,
                    template <typename, typename> class> class M>
hpx_action_t Tree<S, T, E, M>::target_partition_ = HPX_ACTION_NULL;


} // namespace dashmm


#endif // __DASHMM_TREE_H__
