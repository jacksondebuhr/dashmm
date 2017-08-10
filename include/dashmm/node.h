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


#ifndef __DASHMM_NODE_H__
#define __DASHMM_NODE_H__


/// \file
/// \brief tree Node related types


// C library
#include <cstdlib>
#include <cmath>
#include <cstdint>
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
#include "dashmm/index.h"


namespace dashmm {


// Forward declare Registrars for the objects in this file
template <typename Record>
class NodeRegistrar;


/// A Node of a tree.
///
/// This is a template over the record type. In DASHMM, this will be either
/// the Source type or the Target type.
template <typename Record>
class Node {
 public:
  using record_t = Record;
  using node_t = Node<Record>;
  using arrayref_t = ArrayRef<Record>;

  /// The default constructor allocates nothing, and sets all to zero.
  Node() : idx{}, parts{}, parent{nullptr}, dag{this}, first_{0} {
    for (int i = 0; i < 8; ++i) {
      child[i] = nullptr;
    }
    sema_ = HPX_NULL;
    complete_ = HPX_NULL;
  }

  /// Constuct with a known index. This will set the index of the node,
  /// and will allocate the completion detection LCO.
  ///
  /// \param index - the index of the node
  Node(Index index)
      : idx{index}, parts{}, parent{nullptr}, dag{index}, first_{0} {
    for (int i = 0; i < 8; ++i) {
      child[i] = nullptr;
    }
    sema_ = HPX_NULL;
    complete_ = hpx_lco_and_new(8);
  }

  /// Construct with an index, a particle segment and a parent
  ///
  /// This will also create the completion detection LCO.
  ///
  /// \param index - the node index
  /// \param parts - the particles inside the volume represented by this node
  /// \param parent - the parent of this node
  Node(Index index, arrayref_t parts, node_t *parent)
      : idx{index}, parts{parts}, parent{parent}, dag{this, index}, first_{0} {
    for (int i = 0; i < 8; ++i) {
      child[i] = nullptr;
    }
    sema_ = HPX_NULL;
    complete_ = hpx_lco_and_new(8);
  }

  /// Any allocated LCOs are cleaned up explicitly.
  ///
  /// The completion LCOs are cleaned up as they finish their detection, and
  /// pass on any messages dependent on them. Very few nodes of the tree ever
  /// get a lock, so those are cleaned up when needed.
  ~Node() { }

  /// Returns the first record into which new records may be copied.
  ///
  /// NOTE: the node should be locked before using this to guarantee consistent
  /// access.
  size_t first() const {return first_;}

  /// Returns the LCO signaling completion of partitioning for this node
  hpx_addr_t complete() const {return complete_;}

  /// Gives the number of particles in the segment owned by this node
  size_t num_parts() const {return parts.n();}

  /// Create a completion detection and gate
  void add_completion() {
    assert(complete_ == HPX_NULL);
    complete_ = hpx_lco_and_new(8);
    assert(complete_ != HPX_NULL);
  }

  /// Create semaphore as needed
  void add_lock() {
    assert(sema_ == HPX_NULL);
    sema_ = hpx_lco_sema_new(1);
    assert(sema_ != HPX_NULL);
  }

  /// Delete semaphore
  void delete_lock() {
    assert(sema_ != HPX_NULL);
    hpx_lco_delete_sync(sema_);
    sema_ = HPX_NULL;
  }

  /// Lock the semaphore
  void lock() const {
    assert(sema_ != HPX_NULL);
    hpx_lco_sema_p(sema_);
  }

  /// Unlock the semaphore
  void unlock() const {
    assert(sema_ != HPX_NULL);
    hpx_lco_sema_v(sema_, HPX_NULL);
  }

  /// Increment the first record that is open
  ///
  /// See first() above. This will increment the location in which new records
  /// can be added to the rebalanced portion of the point data for this node.
  /// This will return if the given increment fills up the remaining spots in
  /// the segment of global array.
  ///
  /// NOTE: This should only be called once this node has been lock()-ed
  ///
  /// NOTE: This is only used for the nodes of the uniform grid.
  ///
  /// \param incr - the number of records that have just been added to the
  ///               segment.
  ///
  /// \returns - true if the segment is full; false otherwise
  bool increment_first(size_t incr) {
    first_ += incr;
    return first_ >= parts.n();
  }

  /// Partition the node
  ///
  /// Perform a partition action with the given synchronization LCO. This can
  /// be HPX_NULL, in which case this will be an asynchronous action invocation
  /// of the partition action.
  ///
  /// \param sync - synchronization LCO to trigger once action is complete
  /// \param threshold - partitioning threshold
  /// \param px - x position of low corner of domain
  /// \param py - y position of low corner of domain
  /// \param pz - z position of low corner of domain
  /// \param size - size of domain
  /// \param same_sandt - nonzero if this is a case where S==T
  void partition(hpx_addr_t sync, int threshold,
                 double px, double py, double pz, double size,
                 int same_sandt) {
    node_t *thisarg = this;
    hpx_call(HPX_HERE, partition_node_, sync,
             &thisarg, &px, &py, &pz, &size, &threshold, &same_sandt);
  }

  /// Partition the node when the given LCO triggers
  ///
  /// Perform a partition action with the given synchronization LCO. This can
  /// be HPX_NULL, in which case this will be an asynchronous action invocation
  /// of the partition action.
  ///
  /// \param when - gate LCO on which to depend before starting the action
  /// \param sync - synchronization LCO to trigger once action is complete
  /// \param threshold - partitioning threshold
  /// \param px - x position of low corner of domain
  /// \param py - y position of low corner of domain
  /// \param pz - z position of low corner of domain
  /// \param size - size of domain
  /// \param same_sandt - nonzero if this is a case where S==T
  void partitionWhen(hpx_addr_t when, hpx_addr_t sync, int threshold,
                     double px, double py, double pz, double size,
                     int same_sandt) {
    node_t *thisarg = this;
    hpx_call_when(when, partition_node_, sync,
                  &thisarg, &px, &py, &pz, &size, &threshold, &same_sandt);
  }

  /// Return the size of the branch below this node
  ///
  /// This will return the total number of descendants of this node.
  ///
  /// \returns - the number of nodes in the branch below this node
  int n_descendants() const {
    // NOTE: the implementation is not recursive because HPX-5 can work with
    // small stack sizes.
    std::vector<const node_t *> V{this};
    int count = 0;
    while (!V.empty()) {
      count += V.size();
      std::vector<const node_t *> C{};
      for (size_t i = 0; i < V.size(); ++i) {
        for (int j = 0; j < 8; ++j) {
          if (V[i]->child[j]) {
            C.push_back(V[i]->child[j]);
          }
        }
      }
      V = std::move(C);
    }
    return count;
  }

  /// Return the number of immediate descendants of this node
  ///
  /// \returns - the number of non-null children of this node.
  int n_children() const {
    int retval{0};
    for (int i = 0; i < 8; ++i) {
      if (child[i] != nullptr) {
        ++retval;
      }
    }
    return retval;
  }

  /// Predicate for testing if this node is a leaf
  ///
  /// \returns - true if a leaf; false otherwise
  bool is_leaf() const {return n_children() == 0;}

  /// Compress the branch information into the provided buffers
  ///
  /// This will compress the branch information, which can then be sent
  /// to other ranks, where is can be extracted into the needed nodes.
  ///
  /// The structure of these is as two arrays, with one entry for each node
  /// except the root of the branch (which is the uniform grid node below
  /// which this branch is found, which is also this). One array gives the
  /// index of the parent of the current node, and the other gives which child
  /// of that parent is the current node.
  ///
  /// \param branch - which child of the parent indicated by tree is this node
  /// \param tree - which node is the parent of the given node
  void compress(int *branch, int *tree) const {
    // NOTE: the non-recursive implementation is to take into account HPX-5's
    // small default stack size.
    std::vector<const node_t *> V{this};
    std::vector<int> V_idx{-1};
    int curr = 0;
    while (!V.empty()) {
      std::vector<const node_t *> C{};
      std::vector<int> C_idx{};
      for (size_t i = 0; i < V.size(); ++i) {
        assert(V[i] != nullptr);
        for (int j = 0; j < 8; ++j) {
          if (V[i]->child[j]) {
            branch[curr] = j;
            tree[curr] = V_idx[i];
            C.push_back(V[i]->child[j]);
            C_idx.push_back(curr++);
          }
        }
      }
      V = std::move(C);
      V_idx = std::move(C_idx);
    }
  }

  /// Extract the branch information from the provided buffers
  ///
  /// This will extract the branch information received from other ranks,
  /// which can then be reified into the needed nodes.
  ///
  /// \param branch - buffers holding the tree structure
  /// \param tree -
  /// \param n_nodes - the number of nodes
  void extract(const int *branch, const int *tree, int n_nodes) {
    // Extract a compressed remote tree representation. As the tree is remote,
    // only {parent, child, idx} fields are needed.

    node_t *descendants = new node_t[n_nodes]{};

    // The compressed tree is created in depth first fashion. And there are two
    // choices here to fill in the parent, child, and idx fields of the
    // descendants.

    // Approach I: Setup parent and child, which is an embarassingly parallel
    // operation on the @p branch and @p tree. Afterwards, fan out along the
    // tree to fill in idx.

    // Approach II: Go over the input sequentially. For each node encountered,
    // by the depth first property, the index of its parent is already set.
    // So one can finish in one loop.

    // If on each rank, there are multiple subtrees being merged, approach II
    // might be sufficient. Approach I can be considered if finer granularity is
    // needed.

    // Approach II is implemented here.
    for (int i = 0; i < n_nodes; ++i) {
      int pos = tree[i];
      int which = branch[i];
      node_t *curr = &descendants[i];
      node_t *parent = (pos < 0 ? this : &descendants[pos]);

      curr->parent = parent;
      curr->idx = parent->idx.child(which);
      curr->dag.set_index(parent->idx.child(which));
      parent->child[which] = curr;
    }
  }

  /// Remove a given child from this node
  void remove_child(const node_t *target) {
    bool found{false};
    for (int i = 0; i < 8; ++i) {
      if (child[i] == target) {
        child[i] = nullptr;
        found = true;
        break;
      }
    }
    assert(found);
  }

  /// This will recurse down to the uniform level and remove pointless nodes
  ///
  /// This returns true if this node has no children after the work is
  /// completed for this node.
  ///
  /// NOTE: This should be a safe recursion because the depth of the top
  /// part of the tree grows very slowly.
  bool removeDownwardLinks(int limit, int level) {
    if (level < limit) {
      for (int i = 0; i < 8; ++i) {
        if (child[i] && child[i]->removeDownwardLinks(limit, level + 1)) {
          child[i]->parent = nullptr; // Remove upward link
          child[i] = nullptr;         // Remove downward link
        }
      }
    }
    return (is_leaf() && parts.n() == 0);
  }

  /// Destroy the data in a branch
  ///
  /// This is only ever called on the local nodes. This will traverse the
  /// branch and delete all nodes below the given node.
  ///
  /// \param root - the root of the branch.
  static void destroy_branch(node_t *root) {
    // NOTE: this is not recursive because HPX-5 has small default stacks.
    std::vector<node_t *> V{root};
    while (!V.empty()) {
      std::vector<node_t *> C{};
      for (size_t i = 0; i < V.size(); ++i) {
        assert(V[i] != nullptr);
        for (int j = 0; j < 8; ++j) {
          if (V[i]->child[j]) {
            C.push_back(V[i]->child[j]);
          }
        }
        delete V[i];
        V[i] = nullptr;
      }
      V = std::move(C);
    }
  }


  // TODO: These should likely all be privatized.
  Index idx;                      /// index of the node
  arrayref_t parts;               /// segment for this node
  node_t *parent;                 /// parent node
  node_t *child[8];               /// children of this node
  DAGInfo dag;                    /// The DAG info for this node


 private:
  friend class NodeRegistrar<Record>;

  /// Partition the node
  ///
  /// In addition to sorting the points associated with this node, this will
  /// create the needed children and schedule the work of partitioning for
  /// those children.
  ///
  /// @p same_sandt will be nonzero only for target nodes, and only sometimes.
  /// In this case, the following does not sort the records, but merely finds
  /// the split point, which will have been established already by the source
  /// tree partioning.
  ///
  /// \param n - the tree node in question
  /// \param px - x position of low corner of domain
  /// \param py - y position of low corner of domain
  /// \param pz - z position of low corner of domain
  /// \param size - size of domain
  /// \param threshold - the partitioning threshold
  /// \param same_sandt - is this a run where the sources and targets are
  ///                     identical.
  static int partition_node_handler(node_t *n,
                                    double px, double py, double pz,
                                    double size,
                                    int threshold,
                                    int same_sandt) {
    DomainGeometry geo{Point{px, py, pz}, size};

    size_t num_points = n->num_parts();
    assert(num_points >= 1);
    bool is_leaf = num_points <= (size_t)threshold;

    if (n->parent && n->parent->complete() != HPX_NULL) {
      hpx_call_when_with_continuation(n->complete_,
          n->parent->complete(), hpx_lco_set_action,
          complete_, hpx_lco_delete_action,
          nullptr, 0);
    }

    if (is_leaf) {
      // No children means this node is done with partitioning
      hpx_lco_and_set_num(n->complete_, 8, HPX_NULL);
    } else {
      // Compute center of the node
      double h = geo.size() / pow(2, n->idx.level());
      double center_x = geo.low().x() + (n->idx.x() + 0.5) * h;
      double center_y = geo.low().y() + (n->idx.y() + 0.5) * h;
      double center_z = geo.low().z() + (n->idx.z() + 0.5) * h;

      // Get the local data
      record_t *p = n->parts.data();

      // Sort the particles among the children
      record_t *splits[9]{};
      splits[0] = p;
      splits[8] = &p[num_points];

      auto z_comp = [&center_z](record_t &a) {
        return a.position.z() < center_z;
      };
      auto y_comp = [&center_y](record_t &a) {
        return a.position.y() < center_y;
      };
      auto x_comp = [&center_x](record_t &a) {
        return a.position.x() < center_x;
      };

      if (same_sandt) {
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

      // Perform some counting
      int stat[8]{};
      int offset[8]{};

      offset[0] = 0;
      stat[0] = splits[1] - splits[0];
      for (int i = 1; i < 8; ++i) {
        stat[i] = splits[i + 1] - splits[i];
        offset[i] = offset[i - 1] + stat[i - 1];
      }

      // Create child nodes
      for (int i = 0; i < 8; ++i) {
        if (stat[i]) {
          auto cparts = n->parts.slice(offset[i], stat[i]);
          node_t *cnd = new node_t{n->idx.child(i), cparts, this};
          n->child[i] = cnd;
          cnd->partition(HPX_NULL, threshold, px, py, pz, size, same_sandt);
        } else {
          hpx_lco_and_set(n->complete_, HPX_NULL);
        }
      }
    }

    return HPX_SUCCESS;
  }

  size_t first_;            /// first record that is available
  hpx_addr_t sema_;         /// restrict concurrent modification
  hpx_addr_t complete_;     /// This is used to indicate that partitioning is
                            ///  complete

  static hpx_action_t partition_node_;
};

template <typename R>
hpx_action_t Node<R>::partition_node_ = HPX_ACTION_NULL;


} // dashmm


#endif // __DASHMM_NODE_H__
