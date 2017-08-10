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


#ifndef __DASHMM_TREE_H__
#define __DASHMM_TREE_H__


/// \file
/// \brief Tree related types


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


namespace dashmm {


// Forward declare registrar for Tree object
template <typename Source, typename Target, typename Record>
class TreeRegistrar;


/// The tree class represents one of the trees in a given evaluation
///
/// The Tree is a template over the typical DASHMM types, as well as a
/// Record type that indicates if this is a Source tree or a Target tree.
///
/// The nodes are arranged in a hybrid way in this tree. The top of the tree
/// (up to and including the finest uniform level) are allocated in an array.
/// The branches owned by this rank are allocated in the traditional fashion,
/// one node at a time. Branches from remotes are allocated in an array.
/// The two different schemes reflect that we sometimes know how many nodes we
/// shall need, and other times we do not.
template <typename Source, typename Target, typename Record>
class Tree {
 public:
  using record_t = Record;
  using node_t = Node<Record>;
  using sourcenode_t = Node<Source>;
  using targetnode_t = Node<Target>;
  using arrayref_t = ArrayRef<Record>;
  using tree_t = Tree<Source, Target, Record>;
  using sourcetree_t = Tree<Source, Target, Source>;

  /// Tree construction just default initializes the object
  Tree() : root_{nullptr}, unif_grid_{nullptr}, unif_done_{HPX_NULL},
           sorted_{} { }

  arrayref_t sorted() const {return sorted_;}

  // TODO: I do not like this. See dualtree_t::create_DAG for the use case.
  // also dualtree_t::collect_DAG_nodes.
  // And dualtree_t::create_expansions_from_DAG
  node_t *root() {return root_;}

  /// Setup some basic information during initial tree construction
  ///
  /// This action is the target of a broadcast. The basic information about
  /// the tree is set up in this action. The most important of which is the
  /// allocation of the top portions of the tree containing the root, and
  /// extending to the finest uniform level of partitioning.
  ///
  /// \param sync - the address of an LCO to use for synchronization
  /// \param unif_level - the uniform partitioning level
  void setupBasics(hpx_addr_t sync, int unif_level) {
    tree_t *thisarg = this;
    hpx_call(HPX_HERE, setup_basics_, sync, &thisarg, &unif_level);
  }

  /// Delete the data stored in this tree
  ///
  /// This is an asynchronous operation. It will return after the parecel is
  /// sent. To track completion of the action, provided a non-HPX_NULL @p sync
  /// argument to this method.
  ///
  /// \param sync - the synchronization LCO
  /// \param global_tree - the global Tree
  /// \param ndim - the number of uniform level nodes
  /// \param rmap - the mapping from those nodes to owning rank
  static void delete_tree(hpx_addr_t sync,
                          RankWise<tree_t> global_tree,
                          int ndim,
                          int *rmap) {
    hpx_addr_t targ = global_tree.data();
    hpx_call(HPX_HERE, delete_tree_, sync, &targ, &ndim, &rmap);
  }

  /// Assign points to the uniform grid.
  ///
  /// This will assign the points to the uniform grid. This gives the points
  /// the id (in the Morton Key sense) of the box to which they are assigned,
  /// and it will count the numbers in each box. This is a key first step to
  /// computing the distribution of the sources and targets during tree
  /// construction.
  ///
  /// \param P - the records
  /// \param npts - the number of records
  /// \param geo - the overall domain geometry
  /// \param unif_level - the uniform partitioning level
  /// \param gid [out] - the morton key for each record
  /// \param count [out] - the number of points per uniform grid node
  static void assign_points_to_unif_grid(hpx_addr_t sync,
                                         const record_t *P,
                                         int npts,
                                         const DomainGeometry *geo,
                                         int unif_level,
                                         int *gid,
                                         int *count) {
    hpx_call(HPX_HERE, assign_points_, sync,
             &P, &npts, &geo, &unif_level, &gid, &count);
  }

  /// Reorder the particles according to their place in the uniform grid
  ///
  /// This will rearrange the particles into their bin order.
  ///
  /// \param p_in - the input records; will be sorted
  /// \param npts - the number of records
  /// \param dim3 - the size of the uniform grid
  /// \param gid_of_points - the morton key for the records
  /// \param count - the number of records per bin
  /// \param retval [out] - offsets into the record list for each
  static void group_points_on_unif_grid(hpx_addr_t sync,
                                        record_t *p_in,
                                        int npts,
                                        int dim3,
                                        int *gid_of_points,
                                        const int *count,
                                        int **retval) {
    hpx_call(HPX_HERE, group_points_, sync,
             &p_in, &npts, &dim3, &gid_of_points, &count, &retval);
  }

  /// Send a compressed version of the node to the other ranks
  ///
  /// This will send a compressed version of the branch of the tree below
  /// the given uniform node to the other ranks. This is done to make it
  /// so that each rank can have a complete view of the structure of the
  /// tree.
  ///
  /// \param sync - the LCO to set once the action is done
  /// \param n_idx - the index of the uniform level node
  /// \param rwtree - the global address of the Tree
  void sendNode(hpx_addr_t sync, int n_idx, hpx_addr_t rwtree) {
    node_t *curr = &unif_grid_[n_idx];
    hpx_call_when(curr->complete(), HPX_HERE, send_node_, sync,
                  &curr, &n_idx, &rwtree);
  }

  /// Initialize point exchange between ranks
  ///
  /// This is largely just a setup procedure to compute a few things that are
  /// needed for point exchanges. The most important of which is the global
  /// offset in this rank's data where each incoming batch of points will be
  /// stored in the sorted array. Also, this allocates the segment of memory
  /// that will store the sorted data.
  ///
  /// \param dim3 - the number of uniform level nodes
  /// \param rank_map - the node->rank map for the uniform nodes
  /// \param global_count - the global counts in each node of the uniform grid
  /// \param local_count - the local counts in each node of the uniform grid
  /// \param local_offset - where in the local data are the points for each node
  ///                       of the uniform grid
  /// \param geo - the domain geometry for the tree
  /// \param threshold - the partitioning threshold for the tree
  /// \param temp - the local point data
  void initPointExchange(int dim3,
                         int *rank_map,
                         const int *global_count,
                         const int *local_count,
                         const int *local_offset,
                         const DomainGeometry &geo,
                         int threshold,
                         const record_t *temp) {
    arrayref_t sorted_ref{};

    // Compute global_offset
    int my_rank = hpx_get_my_rank();
    int *global_offset = new int[dim3]();
    size_t num_points = rank_map[0] == my_rank ? global_count[0] : 0;
    for (int i = 1; i < dim3; ++i) {
      num_points += rank_map[i] == my_rank ? global_count[i] : 0;
      int increment = rank_map[i - 1] == my_rank ? global_count[i - 1] : 0;
      global_offset[i] = global_offset[i - 1] + increment;
    }

    if (num_points > 0) {
      record_t *sorted_data = new record_t[num_points];
      sorted_ref = arrayref_t{sorted_data, num_points};

      for (int i = 0; i < dim3; ++i) {
        if (rank_map[i] != my_rank) continue;

        node_t *curr = &unif_grid_[i];
        curr->parts = sorted_ref.slice(global_offset[i], global_count[i]);

        if (local_count[i]) {
          // Copy local points before merging remote points
          curr->lock();
          auto sorted = curr->parts.data();
          for (size_t j = 0; j < local_count[i]; ++j) {
            sorted[j + curr->first()] = temp[local_offset[i] + j];
          }

          if (curr->increment_first(local_count[i])) {
            // This grid does not expect remote points.
            // Spawn adaptive partitioning
            Point geo_pt = geo.low();
            double px = geo_pt.x();
            double py = geo_pt.y();
            double pz = geo_pt.z();
            double size = geo.size();
            curr->partition(HPX_NULL, threshold, px, py, pz, size, 0);
          }
          curr->unlock();
        }
      }
    }
    delete [] global_offset;

    sorted_ = sorted_ref;
    hpx_lco_and_set(unif_done_, HPX_NULL);
  }

  /// Initialize point exchange between ranks for S == T case
  ///
  /// This is largely just a setup procedure to compute a few things that are
  /// needed for point exchanges. The most important of which is the global
  /// offset in this rank's data where each incoming batch of points will be
  /// stored in the sorted array. Also, this allocates the segment of memory
  /// that will store the sorted data.
  ///
  /// \param tree - the tree object
  /// \param dim3 - the number of uniform nodes
  /// \param rank_map - the node->rank map for the uniform nodes
  /// \param global_count - the global counts in each node of the uniform grid
  /// \param local_count - the local counts in each node of the uniform grid
  /// \param local_offset - where in the local data are the points for each node
  ///                       of the uniform grid
  /// \param geo - the domain geometry for the tree
  /// \param threshold - the partitioning threshold for the tree
  /// \param temp - the local point data
  /// \param n - the uniform grid nodes
  /// \param source_tree - the source tree
  void initPointExchangeSameSAndT(int dim3,
                                  int *rank_map,
                                  const int *global_count,
                                  const int *local_count,
                                  const int *local_offset,
                                  const DomainGeometry *geo,
                                  int threshold,
                                  const record_t *temp,
                                  sourcetree_t *source_tree) {
    int my_rank = hpx_get_my_rank();
    size_t num_points{0};
    for (int i = 0; i < dim3; ++i) {
      if (rank_map[i] == my_rank) {
        num_points += global_count[i];
      }
    }

    if (num_points > 0) {
      for (int i = 0; i < dim3; ++i) {
        if (rank_map[i] != my_rank) continue;

        // When S == T, we need to set the parts on the target to be the
        // parts on the equivalent source node
        node_t *curr = &unif_grid_[i];
        sourcenode_t *curr_source = &source_tree->unif_grid_[i];
        ArrayRef<Source> sparts = curr_source->parts;
        curr->parts = arrayref_t{(record_t *)sparts.data(), sparts.n()};

        if (local_count[i]) {
          curr->lock();
          // still need to increment first, and still need to partition
          // if things are ready. The only catch is that we have to wait for
          // the equivalent source node to be ready before we start.
          if (curr->increment_first(local_count[i])) {
            // This grid does not expect remote points.
            // Spawn adaptive partitioning
            Point geo_pt = geo->low();
            double px = geo_pt.x();
            double py = geo_pt.y();
            double pz = geo_pt.z();
            double size = geo->size();
            curr->partitionWhen(curr_source->complete(), HPX_NULL,
                                threshold, px, py, pz, size, 1);
          }
          curr->unlock();
        }
      }
    }

    // Again, the target tree reuses the data from the source tree
    ArrayRef<Source> ssort = source_tree->sorted();
    sorted_ = arrayref_t{(record_t *)ssort.data(), ssort.n()};
    hpx_lco_and_set(unif_done_, HPX_NULL);
  }

  /// Merge incoming points into the local array
  ///
  /// This action merges particular points with the sorted list. Also, if this
  /// is the last set of points that are merged, this will go ahead and start
  /// the adaptive partitioning of that part of the local tree.
  ///
  /// This is an asynchronous call. To track completion of the actino launched
  /// by this method, provide an LCO that will be set as @p sync. If no
  /// such tracking is needed, the caller may provide HPX_NULL.
  ///
  /// \param sync - LCO to trigger at the completion of the action
  /// \param temp - the local data
  /// \param n_idx - the uniform grid node's index
  /// \param n_arrived - the number arriving in this message
  /// \param geo - the domain geometry
  /// \param thresh - the partitioning threshold
  void mergePoints(hpx_addr_t sync,
                   record_t *temp,
                   int n_idx,
                   int n_arrived,
                   const DomainGeometry *geo,
                   int thresh) {
    node_t *node = &unif_grid_[n_idx];
    Point geo_pt = geo->low();
    double px = geo_pt.x();
    double py = geo_pt.y();
    double pz = geo_pt.z();
    double size = geo->size();
    hpx_call(HPX_HERE, merge_points_, sync,
             &temp, &node, &n_arrived, &px, &py, &pz, &size, &thresh);
  }

  /// Merge incoming points into the local array
  ///
  /// This performs the S == T version of point merging. In this version,
  /// no points are actually merged, as that is handles in the Source version of
  /// this routine. Instead, this merely updates the first counter, and
  /// calls to partitioning when appropriate. Note, this is now a call when
  /// waiting on the completion of the partitioning in the source tree.
  ///
  /// This spawns the action asynchronously, and completion of the action can
  /// be detected with the @p sync LCO argument. If no such completeion
  /// detection is needed, the caller may provide HPX_NULL.
  ///
  /// \param sync - completion detection LCO
  /// \param n_idx - the index of the uniform node in question
  /// \param n_arrived - the number arriving in this message
  /// \param stree - the source tree
  /// \param geo - the domain geometry
  /// \param thresh - the partitioning threshold
  void mergePointsSameSAndT(hpx_addr_t sync,
                            int n_idx,
                            int n_arrived,
                            sourcetree_t *stree,
                            const DomainGeometry *geo,
                            int thresh) {
    targetnode_t *tnode = &unif_grid_[n_idx];
    sourcenode_t *snode = &stree->unif_grid_[n_idx];
    Point geo_pt = geo->low();
    double px = geo_pt.x();
    double py = geo_pt.y();
    double pz = geo_pt.z();
    double size = geo->size();
    hpx_call(HPX_HERE, merge_points_same_s_and_t_, sync,
             &tnode, &n_arrived, &snode, &px, &py, &pz, &size, &thresh);
  }

  /// Split the bits of an integer to be used in a Morton Key
  static uint64_t split(unsigned k) {
    uint64_t split = k & 0x1fffff;
    split = (split | split << 32) & 0x1f00000000ffff;
    split = (split | split << 16) & 0x1f0000ff0000ff;
    split = (split | split << 8)  & 0x100f00f00f00f00f;
    split = (split | split << 4)  & 0x10c30c30c30c30c3;
    split = (split | split << 2)  & 0x1249249249249249;
    return split;
  }

  /// Compute the Morton key for a gives set of indices
  static uint64_t morton_key(unsigned x, unsigned y, unsigned z) {
    uint64_t key = 0;
    key |= split(x) | split(y) << 1 | split(z) << 2;
    return key;
  }

  /// Compute the simple 3d key of an index
  ///
  /// \param x - the x component
  /// \param y - the y component
  /// \parma z - the z component
  /// \param max - the max value of any component at the uniform refinement
  ///              level. This will be 2^unif_level
  static uint64_t simple_key(unsigned x, unsigned y, unsigned z, unsigned max) {
    return x + y * max + z * max * max;
  }

  /// Compute the uniform grid index for a given Index
  ///
  /// \param idx - the index of the node
  /// \param uniflevel - the uniform level of the tree
  ///
  /// \returns - the index in the uniform grid that gives the ancestor of
  ///            the given Index
  static int get_unif_grid_index(const Index &idx, int uniflevel) {
    int delta = idx.level() - uniflevel;

    if (delta < 0) {
      return -1;
    }

    int max = 1 << uniflevel;

    if (delta > 0) {
      Index tester = idx.parent(delta);
      return simple_key(tester.x(), tester.y(), tester.z(), max);
    } else {
      return simple_key(idx.x(), idx.y(), idx.z(), max);
    }
  }

  /// Find the LCO address for a given index and a given operation
  ///
  /// \param idx - the Index of the node in question
  /// \param op - the edge type connecting to the index in question
  ///
  /// \returns - global address of the LCO serving as target of the edge
  hpx_addr_t lookup_lco_addx(Index idx, Operation op) {
    // This should walk to the node containing the LCO we care about
    node_t *curr = root_;
    bool not_found = true;
    while (not_found) {
      int dlevel = idx.level() - curr->idx.level();
      assert(dlevel >= 0);
      if (dlevel == 0) {
        break;
      }
      assert(curr->idx == idx.parent(dlevel));
      int which = idx.parent(dlevel - 1).which_child();
      curr = curr->child[which];
    }
    assert(curr->idx == idx);

    hpx_addr_t retval{HPX_NULL};
    switch (op) {
      case Operation::Nop:
        assert(0 && "problem in lookup");
        break;
      case Operation::StoM:  // NOTE: fallthrough here
      case Operation::StoL:
      case Operation::MtoM:
      case Operation::MtoL:
      case Operation::LtoL:
        assert(curr->dag.has_normal());
        retval = curr->dag.normal()->global_addx;
        break;
      case Operation::MtoT: // NOTE: fallthrough here
      case Operation::LtoT:
      case Operation::StoT:
        assert(curr->dag.has_parts());
        retval = curr->dag.parts()->global_addx;
        break;
      case Operation::MtoI: // NOTE: fallthrough
      case Operation::ItoI:
        assert(curr->dag.has_interm());
        retval = curr->dag.interm()->global_addx;
        break;
      case Operation::ItoL:
        assert(curr->dag.has_normal());
        retval = curr->dag.normal()->global_addx;
        break;
    }

    return retval;
  }

  /// Wait for the uniform level work to be finished
  ///
  /// This call will block the calling HPX-5 thread.
  void waitForUnifDone() {
    hpx_lco_wait(unif_done_);
  }

  /// Returns the uniform refinement level node's completion LCO
  hpx_addr_t unifNodeCompletion(int i) {
    return unif_grid_[i].complete();
  }

  /// Prune unused topnodes from the tree
  ///
  /// Topnodes are those at the uniform refinement level or lower. This will
  /// remove the links from the root to these nodes so that various traversals
  /// will not visit nodes with no sources or targets, and will thus save
  /// some overhead of HPX-5 actions that do nothing.
  ///
  /// \param counts - the counts per node of the unif level nodes
  /// \param dim3 - the number of uniform level nodes
  /// \param unif_level - the uniform level in the tree
  void pruneTopnodes(const int *counts, int dim3, int unif_level) {
    for (int i = 0; i < dim3; ++i) {
      node_t *curr = &unif_grid_[i];
      if (counts[i] == 0) {
        curr->parent->remove_child(curr);
      }
    }

    root_->removeDownwardLinks(unif_level - 1, 0);
  }

private:
  // NOTE: One of these is superfluous
  friend class Tree<Source, Target, Source>;
  friend class Tree<Source, Target, Target>;

  friend class TreeRegistrar<Source, Target, Record>;

  /// Setup some basic information during initial tree construction
  ///
  /// This action is the target of a broadcast. The basic information about
  /// the tree is set up in this action. The most important of which is the
  /// allocation of the top portions of the tree containing the root, and
  /// extending to the finest uniform level of partitioning.
  ///
  /// \param tree - the address of the tree on which to act
  /// \param unif_level - the uniform partitioning level
  static int setup_basics_handler(tree_t *tree, int unif_level) {
    tree->unif_done_ = hpx_lco_and_new(1);
    assert(tree->unif_done_ != HPX_NULL);

    // Setup unif_grid
    int n_top_nodes{1};
    for (int i = 1; i < unif_level; ++i) {
      n_top_nodes += pow(8, i);
    }
    int dim3 = pow(8, unif_level);
    tree->root_ = new node_t[n_top_nodes + dim3]{};
    tree->unif_grid_ = &tree->root_[n_top_nodes];

    tree->root_[0].idx = Index{0, 0, 0, 0};
    int startingnode{0};
    int stoppingnode{1};
    for (int level = 0; level < unif_level ; ++level) {
      for (int nd = startingnode; nd < stoppingnode; ++nd) {
        node_t *snode = &tree->root_[nd];

        if (level != unif_level - 1) {
          // At the lower levels, we do not require the ordering
          int firstchild = nd * 8 + 1;
          for (int i = 0; i < 8; ++i) {
            node_t *scnode = &tree->root_[firstchild + i];
            snode->child[i] = scnode;
            scnode->idx = snode->idx.child(i);
            scnode->dag.set_index(snode->idx.child(i));
            scnode->parent = snode;
          }
        } else {
          // The final level does need the ordering
          for (int i = 0; i < 8; ++i) {
            Index cindex = snode->idx.child(i);
            uint64_t kval = simple_key(cindex.x(),
                                       cindex.y(),
                                       cindex.z(),
                                       1 << unif_level);
            snode->child[i] = &tree->unif_grid_[kval];
            tree->unif_grid_[kval].idx = cindex;
            tree->unif_grid_[kval].dag.set_index(cindex);
            tree->unif_grid_[kval].parent = snode;
            tree->unif_grid_[kval].add_lock();
            tree->unif_grid_[kval].add_completion();
          }
        }
      }

      startingnode += pow(8, level);
      stoppingnode += pow(8, level + 1);
    }

    return HPX_SUCCESS;
  }

  /// Assign points to the uniform grid.
  ///
  /// This will assign the points to the uniform grid. This gives the points
  /// the id (in the Morton Key sense) of the box to which they are assigned,
  /// and it will count the numbers in each box. This is a key first step to
  /// computing the distribution of the sources and targets during tree
  /// construction.
  ///
  /// \param P - the records
  /// \param npts - the number of records
  /// \param geo - the overall domain geometry
  /// \param unif_level - the uniform partitioning level
  /// \param gid [out] - the morton key for each record
  /// \param count [out] - the number of points per uniform grid node
  ///
  /// \return HPX_SUCCESS
  static int assign_points_to_unif_grid_handler(const record_t *P, int npts,
                                                const DomainGeometry *geo,
                                                int unif_level, int *gid,
                                                int *count) {
    Point corner = geo->low();
    double scale = 1.0 / geo->size();

    int dim = pow(2, unif_level);
    for (int i = 0; i < npts; ++i) {
      const record_t *p = &P[i];
      int xid = std::min(dim - 1,
                         (int)(dim * (p->position.x() - corner.x()) * scale));
      int yid = std::min(dim - 1,
                         (int)(dim * (p->position.y() - corner.y()) * scale));
      int zid = std::min(dim - 1,
                         (int)(dim * (p->position.z() - corner.z()) * scale));
      gid[i] = simple_key(xid, yid, zid, dim);
      count[gid[i]]++;
    }

    return HPX_SUCCESS;
  }

  /// Reorder the particles according to their place in the uniform grid
  ///
  /// This will rearrange the particles into their bin order.
  ///
  /// \param p_in - the input records; will be sorted
  /// \param npts - the number of records
  /// \param dim3 - the size of the uniform grid
  /// \param gid_of_points - the morton key for the records
  /// \param count - the number of records per bin
  /// \param retval [out] - offsets into the record list for each
  ///
  /// \returns - HPX_SUCCESS
  static int group_points_on_unif_grid_handler(record_t *p_in,
                                               int npts,
                                               int dim3,
                                               int *gid_of_points,
                                               const int *count,
                                               int **retval) {
    int *offset = new int[dim3]();

    offset[0] = 0;
    for (int i = 1; i < dim3; ++i) {
      offset[i] = offset[i - 1] + count[i - 1];
    }

    // Set gid to the final location of the particle; NOTE: this will modify
    // the offset. This is corrected below
    for (int i = 0; i < npts; ++i) {
      int gid = gid_of_points[i];
      gid_of_points[i] = offset[gid];
      offset[gid]++;
    }

    // Declare some loop variables
    record_t source, save;
    int source_sort, save_sort;
    int isource, isave, idest;

    // O(N) in place rearrangement
    for (int i = 0; i < npts; ++i) {
      while (gid_of_points[i] != i) {
        record_t save = p_in[gid_of_points[i]];
        int idx = gid_of_points[gid_of_points[i]];

        p_in[gid_of_points[i]] = p_in[i];
        gid_of_points[gid_of_points[i]] = gid_of_points[i];

        gid_of_points[i] = idx;
        p_in[i] = save;
      }
    }

    // Correct offset
    for (int i = 0; i < dim3; ++i) {
      offset[i] -= count[i];
    }

    *retval = offset;

    return HPX_SUCCESS;
  }

  /// Merge incoming points into the local array
  ///
  /// This action merges particular points with the sorted list. Also, if this
  /// is the last set of points that are merged, this will go ahead and start
  /// the adaptive partitioning of that part of the local tree.
  ///
  /// \param temp - the local data
  /// \param n - the uniform grid node
  /// \param n_arrived - the number arriving in this message
  /// \param px - the low corner of the domain geometry
  /// \param py - the low corner of the domain geometry
  /// \param pz - the low corner of the domain geometry
  /// \param size - the size of the domain geometry
  /// \param thresh - the partitioning threshold
  ///
  /// \returns - HPX_SUCCESS
  static int merge_points_handler(record_t *temp, node_t *n, int n_arrived,
                                  double px, double py, double pz, double size,
                                  int thresh) {
    // Note: all the pointers are local to the calling rank.
    n->lock();
    size_t first = n->first();
    record_t *p = n->parts.data();
    std::copy(temp, temp + n_arrived, p + first);

    if (n->increment_first(n_arrived)) {
      n->partition(HPX_NULL, thresh, px, py, pz, size, 0);
    }
    n->unlock();

    return HPX_SUCCESS;
  }

  /// Merge incoming points into the local array
  ///
  /// This action performs the S == T version of point merging. In this version,
  /// no points are actually merged, as that is handles in the Source version of
  /// this routine. Instead, this merely updates the first counter, and
  /// calls to partitioning when appropriate. Note, this is now a call when
  /// waiting on the completion of the partitioning in the source tree.
  ///
  /// \param target_node - the target node in question
  /// \param n_arrived - the number arriving in this message
  /// \param source_node - the source node in question
  /// \param px - the low corner of the domain geometry
  /// \param py - the low corner of the domain geometry
  /// \param pz - the low corner of the domain geometry
  /// \param size - the size of the domain geometry
  /// \param thresh - the partitioning threshold
  ///
  /// \returns - HPX_SUCCESS
  static int merge_points_same_s_and_t_handler(targetnode_t *target_node,
                                               int n_arrived,
                                               sourcenode_t *source_node,
                                               double px,
                                               double py,
                                               double pz,
                                               double size,
                                               int thresh) {
    target_node->lock();
    if (target_node->increment_first(n_arrived)) {
      int ssat = 1;
      target_node->partitionWhen(source_node->complete(), HPX_NULL,
                                 thresh, px, py, pz, size, 1);
    }
    target_node->unlock();

    return HPX_SUCCESS;
  }

  /// This action sends compressed node representation to remote localities.
  ///
  /// This action is called once the individual grids are done.
  ///
  /// \param curr - the uniform grid node
  /// \param id - the index of the uniform grid node in question
  /// \param rwaddr - global address of the Tree
  ///
  /// \returns - HPX_SUCCESS
  static int send_node_handler(node_t *curr, int id, hpx_addr_t rwaddr) {
    // Exclude curr as it is already allocated on remote localities
    int n_nodes = curr->n_descendants() - 1;
    size_t msgsize = sizeof(int) * (2 + n_nodes * 2) + sizeof(hpx_addr_t);
    char *message_buffer = new char[msgsize];
    hpx_addr_t *rwdata = reinterpret_cast<hpx_addr_t *>(message_buffer);
    int *compressed_tree = reinterpret_cast<int *>(
                                message_buffer + sizeof(hpx_addr_t));

    *rwdata = rwaddr;
    compressed_tree[0] = id; // where to merge
    compressed_tree[1] = n_nodes; // # of nodes
    if (n_nodes) {
      int *branch = &compressed_tree[2];
      int *tree = &compressed_tree[2 + n_nodes];
      curr->compress(branch, tree);
    }

    int rank = hpx_get_my_rank();
    int num_ranks = hpx_get_num_ranks();
    hpx_addr_t done = hpx_lco_and_new(num_ranks - 1);

    for (int r = 0; r < num_ranks; ++r) {
      if (r != rank) {
        hpx_parcel_t *p = hpx_parcel_acquire(message_buffer, msgsize);
        hpx_parcel_set_target(p, HPX_THERE(r));
        hpx_parcel_set_action(p, recv_node_);
        hpx_parcel_send(p, done);
      }
    }

    // Clear local memory once all parcels are sent
    hpx_lco_wait(done);
    hpx_lco_delete_sync(done);

    delete [] message_buffer;
    return HPX_SUCCESS;
  }

  /// Receive partitioned tree nodes from a remote locality
  ///
  /// This action receives the compressed representation of a remote branch
  /// of the tree, and then expands it into node objects, connectint it to
  /// the correct portion of the tree at this locality. In general, there will
  /// be one such message per uniform grid node not owned by this locality, per
  /// tree.
  ///
  /// \param message_buffer - the message
  /// \param UNUSED - the size of the buffer
  ///
  /// \returns - HPX_SUCCESS
  static int recv_node_handler(char *message_buffer, size_t UNUSED) {
    hpx_addr_t *rwdata = reinterpret_cast<hpx_addr_t *>(message_buffer);
    int *compressed_tree = reinterpret_cast<int *>(
                                message_buffer + sizeof(hpx_addr_t));
    RankWise<tree_t> global_tree{*rwdata};
    auto local_tree = global_tree.here();
    int id = compressed_tree[0];
    int n_nodes = compressed_tree[1];

    node_t *curr = &local_tree->unif_grid_[id];

    const int *branch = &compressed_tree[2];
    const int *tree = &compressed_tree[2 + n_nodes];
    curr->extract(branch, tree, n_nodes);

    hpx_lco_and_set_num(curr->complete(), 8, HPX_NULL);

    return HPX_SUCCESS;
  }

  /// Destroy allocated data for this tree.
  ///
  /// This will delete the local branches of the tree as well as destroying
  /// any locks allocated for the uniform grid, and will free the remote
  /// branches.
  ///
  /// \param tree - the tree on which to act
  /// \param ndim - the size of the uniform grid
  /// \param rmap - the rank map of the nodes
  static int delete_tree_handler(hpx_addr_t tree, int ndim, int *rmap) {
    RankWise<tree_t> global_tree{tree};
    auto local_tree = global_tree.here();

    // NOTE: The difference here is that there are two different allocation
    // schemes for the nodes.

    int rank = hpx_get_my_rank();
    for (int i = 0; i < ndim; ++i) {
      if (rmap[i] == rank) {
        node_t *curr = &local_tree->unif_grid_[i];
        curr->delete_lock();

        for (int j = 0; j < 8; ++j) {
          if (curr->child[j]) {
            node_t::destroy_branch(curr->child[j]);
          }
        }

        hpx_lco_delete_sync(curr->complete());
      } else {
        node_t *curr = &local_tree->unif_grid_[i];
        curr->delete_lock();

        for (int j = 0; j < 8; ++j) {
          node_t *child = curr->child[j];
          if (child) {
            delete [] child;
            break;
          }
        }
      }
    }

    delete [] local_tree->root_;

    hpx_lco_delete_sync(local_tree->unif_done_);

    return HPX_SUCCESS;
  }

  node_t *root_;            /// Root of the tree
  node_t *unif_grid_;       /// The uniform grid
  hpx_addr_t unif_done_;    /// An LCO indicating that the uniform partition is
                            /// complete
  arrayref_t sorted_;       /// A reference to the sorted point data owned by
                            /// this tree.

  static hpx_action_t setup_basics_;
  static hpx_action_t delete_tree_;
  static hpx_action_t recv_node_;
  static hpx_action_t send_node_;
  static hpx_action_t assign_points_;
  static hpx_action_t group_points_;
  static hpx_action_t merge_points_;
  static hpx_action_t merge_points_same_s_and_t_;
};

template <typename S, typename T, typename R>
hpx_action_t Tree<S, T, R>::setup_basics_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R>
hpx_action_t Tree<S, T, R>::delete_tree_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R>
hpx_action_t Tree<S, T, R>::recv_node_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R>
hpx_action_t Tree<S, T, R>::send_node_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R>
hpx_action_t Tree<S, T, R>::assign_points_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R>
hpx_action_t Tree<S, T, R>::group_points_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R>
hpx_action_t Tree<S, T, R>::merge_points_ = HPX_ACTION_NULL;

template <typename S, typename T, typename R>
hpx_action_t Tree<S, T, R>::merge_points_same_s_and_t_ = HPX_ACTION_NULL;


} // namespace dashmm


#endif // __DASHMM_TREE_H__
