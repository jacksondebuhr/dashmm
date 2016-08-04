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


#ifndef __DASHMM_DAG_H__
#define __DASHMM_DAG_H__


/// \file include/dashmm/dag.h
/// \brief Interface for intermediate representation of DAG


#include <string>
#include <vector>

#include <hpx/hpx.h>

#include "dashmm/index.h"


namespace dashmm {


/// Operation codes to indicate the type of edge
enum class Operation {
  Nop,
  StoM,
  StoL,
  MtoM,
  MtoL,
  LtoL,
  MtoT,
  LtoT,
  StoT,
  MtoI,
  ItoI,
  ItoL
};


struct DAGNode;


/// Edge in the explicit representation of the DAG
struct DAGEdge {
  DAGNode *source;          /// Source node of the edge
  DAGNode *target;          /// Target node of the edge
  Operation op;             /// Operation to perform along edge

  DAGEdge() : source{nullptr}, target{nullptr}, op{Operation::Nop} { }
  DAGEdge(DAGNode *start, DAGNode *end, Operation inop)
    : source{start}, target{end}, op{inop} { }
};


/// Node in the explicit representation of the DAG
struct DAGNode {
  std::vector<DAGEdge> out_edges;   /// these are out edges
  std::vector<DAGEdge> in_edges;    /// these are in edges
  Index idx;                        /// index of the containing node

  int locality;                  /// the locality where this will be placed
  hpx_addr_t global_addx;        /// global address of object serving this node
  size_t other_member;           /// this is either n_digits for an expansion or
                                 /// n_targets for a target lco

  DAGNode(Index i)
      : out_edges{}, in_edges{}, idx{i}, locality{-1}, global_addx{HPX_NULL},
        other_member{0} { }
  void add_out_edge(DAGNode *end, Operation op) {
    out_edges.push_back(DAGEdge{this, end, op});
  }
  void add_in_edge(DAGNode *start, Operation op) {
    in_edges.push_back(DAGEdge{start, this, op});
  }
};


/// DAG object
///
/// This is the explicit representation of the DAG for the particular
/// evaluation. This object is constructed on a single locality, and so does
/// not inhabit the global address space. The DAG object is a template from
/// which the particular evaluation will be created. The distribution policy
/// will operate on the DAG object to compute a distribution of the work
/// around the system.
///
/// To avoid cumbersome interfaces, this object provides direct access to its
/// member data, which are four vectors of DAGNodes: one each for the source
/// nodes of the DAG, the target nodes of the DAG, the nodes of the DAG
/// associated with the nodes of the source tree, and the nods of the DAG
/// associated with the nodes of the target tree.
///
/// Typically, DASHMM users implementing a new Method will work with DAGInfo
/// objects rather than the DAG directly.
class DAG {
 public:
  DAG() : source_leaves{}, source_nodes{}, target_nodes{}, target_leaves{} { }

  /// Print the DAG out in JSON format.
  ///
  /// This will open a file with the given name, and write out the DAG
  /// information to that file is JSON format. See the implementation for
  /// details about what is included.
  void toJSON(std::string fname);

  static bool compare_edge_locality(const DAGEdge &a, const DAGEdge &b) {
    return (a.target->locality < b.target->locality);
  }

  static bool operation_to_target(Operation op) {
    return op == Operation::MtoT || op == Operation::LtoT
                                 || op == Operation::StoT;
  }

  std::vector<DAGNode *> source_leaves;
  std::vector<DAGNode *> source_nodes;
  std::vector<DAGNode *> target_nodes;
  std::vector<DAGNode *> target_leaves;
};


/// DAG information relevant for a given tree node
///
/// Each node of the trees will have a DAGInfo member. This will store the
/// relevant information about the DAG for each tree node. This is the object
/// that Methods will use to construct the DAG. Each DAGInfo object will
/// represent up to three nodes of the DAG: one each for sources/targets,
/// the normal expansion, and an optional intermediate expansion.
///
/// A DAGInfo will create the normal expansion node when constructed. The
/// particle node or the intermediate node will be created on an as-needed
/// basis. Source and Target nodes will be created during tree construction.
/// Intermediate nodes should be added by Methods that need them.
///
/// This object uses HPX-5 to manage the concurrent modification of the DAG,
/// so this object cannot be used outside of an HPX-5 thread. DASHMM users will
/// have no reason to create these objects directly; the library will manage
/// the creation of these objects.
class DAGInfo {
 public:
  /// Construct the DAGInfo
  ///
  /// This will add the normal DAGNode, but not the particle or intermediate
  /// nodes.
  ///
  /// This cannot be used outside of an HPX-5 thread.
  DAGInfo(Index idx) : idx_{idx} {
    normal_ = nullptr;
    interm_ = nullptr;
    parts_ = nullptr;
    lock_ = hpx_lco_sema_new(1);
    assert(lock_ != HPX_NULL);
  }

  /// Destroy the DAGInfo
  ///
  /// This will free any resources acquired by this object. This cannot be used
  /// outside of an HPX-5 thread.
  ~DAGInfo() {
    hpx_lco_delete_sync(lock_);
    if (normal_) {
      delete normal_;
      normal_ = nullptr;
    }
    if (interm_ != nullptr) {
      delete interm_;
      interm_ = nullptr;
    }
    if (parts_ != nullptr) {
      delete parts_;
      parts_ = nullptr;
    }
  }

  /// Add the normal node
  ///
  /// This will add a normal DAG node for the tree node owning this object.
  void add_normal() {
    assert(normal_ == nullptr);
    normal_ = new DAGNode{idx_};
    assert(normal_ != nullptr);
  }

  /// Add an intermediate node
  ///
  /// This will add an intermediate DAG node for the tree node associated
  /// with this object.
  void add_interm() {
    assert(interm_ == nullptr);
    interm_ = new DAGNode{idx_};
    assert(interm_ != nullptr);
  }

  /// Add a particle node
  ///
  /// This will add a particle DAG node for the tree node associated with this
  /// object. This will represent either sources in the source tree, or
  /// targets in the target tree.
  void add_parts() {
    assert(parts_ == nullptr);
    parts_ = new DAGNode{idx_};
    assert(parts_ != nullptr);
  }

  DAGInfo(const DAGInfo &other) = delete;
  DAGInfo &operator=(const DAGInfo &other) = delete;
  DAGInfo(const DAGInfo &&other) = delete;
  DAGInfo &operator=(const DAGInfo &&other) = delete;

  /// Does the tree node have a normal DAG node?
  bool has_normal() const {return normal_ != nullptr;}

  /// Does the tree node have an intermediate DAG node?
  bool has_interm() const {return interm_ != nullptr;}

  /// Does the tree node have either a source or target DAG node?
  bool has_parts() const {return parts_ != nullptr;}

  /// Retrieve the normal DAG node
  DAGNode *normal() const {return normal_;}

  /// Retrieve the intermediate DAG node
  DAGNode *interm() const {return interm_;}

  /// Return the source or target DAG node
  DAGNode *parts() const {return parts_;}

  /// Sets the global data for the normal DAG node
  ///
  /// DASHMM eventually uses the DAG to produce LCO objects that perform the
  /// represented by the DAG. This routine sets the global address of the
  /// LCO serving the normal DAG node.
  ///
  /// \param expand - the expansion LCO represented by this object's normal
  ///                 DAG node.
  void set_normal_expansion(const hpx_addr_t addx, const int acc) {
    normal_->global_addx = addx;
    normal_->other_member = acc;
  }

  /// Sets the global data for the intermediate DAG node
  ///
  /// DASHMM eventually uses the DAG to produce LCO objects that perform the
  /// represented by the DAG. This routine sets the global address of the
  /// LCO serving the intermediate DAG node.
  ///
  /// \param expand - the expansion LCO represented by this object's
  ///                 intermediate DAG node.
  void set_interm_expansion(const hpx_addr_t addx, const int acc) {
    if (interm_ != nullptr) {
      interm_->global_addx = addx;
      interm_->other_member = acc;
    }
  }

  /// Sets the data for a target particle DAG node
  ///
  /// DASHMM eventually uses the DAG to produce LCO objects that perform the
  /// represented by the DAG. This routine sets the global address of the
  /// LCO serving the target DAG node.
  ///
  /// \param targs - the target LCO represented by this object's
  ///                particle DAG node.
  void set_targetlco(const hpx_addr_t addx, const int num) {
    assert(parts_ != nullptr);
    if (parts_ != nullptr) {
      parts_->global_addx = addx;
      parts_->other_member = num;
    }
  }

  /// Sets the data for the source particle DAG node
  ///
  /// Source nodes of the DAG are not associated with an LCO; the source data
  /// is ready at the beginning of the evaluation. However, it is useful for
  /// this object to collect the global address of the sources in question.
  /// This routine collects that information.
  ///
  /// \param src - a reference to the source particles represented by the
  ///              particle DAG node of this object.
  void set_sourceref(const hpx_addr_t addx, const int num) {
    if (parts_ != nullptr) {
      parts_->global_addx = addx;
      parts_->other_member = num;
    }
  }

  /// Sets locality on the particle node
  void set_parts_locality(int loc) {
    if (parts_ != nullptr) {
      parts_->locality = loc;
    }
  }

  /// Sets locality on the normal node
  void set_normal_locality(int loc) {
    if (parts_ != nullptr && normal_ != nullptr) {
      normal_->locality = loc;
    }
  }

  /// Create an S->M link in the DAG
  ///
  /// This will connect the source DAG node of the given object to this object's
  /// normal DAG node with an S->M edge.
  ///
  /// \param source - the DAGInfo object containing the source node in question
  void StoM(DAGInfo *source) {
    assert(source->has_parts());
    assert(has_normal());
    link_nodes(source, source->parts_, this, normal_, Operation::StoM);
  }

  /// Create an S->L link in the DAG
  ///
  /// This will connect the source DAG node of the given object to this object's
  /// normal DAG node with an S->L edge.
  ///
  /// \param source - the DAGInfo object containing the source node in question
  void StoL(DAGInfo *source) {
    assert(source->has_parts());
    assert(has_normal()); 
    link_nodes(source, source->parts_, this, normal_, Operation::StoL);
  }

  /// Create an M->M link in the DAG
  ///
  /// This will connect the normal DAG node of the given object to this object's
  /// normal DAG node with an M->M edge.
  ///
  /// \param source - the DAGInfo object containing the normal node in question
  void MtoM(DAGInfo *source) {
    assert(source->has_normal());
    assert(has_normal());
    link_nodes(source, source->normal_, this, normal_, Operation::MtoM);
  }

  /// Create an M->L link in the DAG
  ///
  /// This will connect the normal DAG node of the given object to this object's
  /// normal DAG node with an M->L edge.
  ///
  /// \param source - the DAGInfo object containing the normal node in question
  void MtoL(DAGInfo *source) {
    assert(source->has_normal());
    assert(has_normal());
    link_nodes(source, source->normal_, this, normal_, Operation::MtoL);
  }

  /// Create an L->L link in the DAG
  ///
  /// This will connect the normal DAG node of the given object to this object's
  /// normal DAG node with an L->L edge.
  ///
  /// \param source - the DAGInfo object containing the normal node in question
  void LtoL(DAGInfo *source) {
    assert(source->has_normal());
    assert(has_normal());
    link_nodes(source, source->normal_, this, normal_, Operation::LtoL);
  }

  /// Create an M->T link in the DAG
  ///
  /// This will connect this objects's normal DAG node to the given object's
  /// particle DAG node with an M->T edge.
  ///
  /// \param source - the DAGInfo object containing the particle node in
  ///                 question
  void MtoT(DAGInfo *target) {
    assert(has_normal());
    assert(target->has_parts());
    link_nodes(this, normal_, target, target->parts_, Operation::MtoT);
  }

  /// Create an L->T link in the DAG
  ///
  /// This will connect this objects's normal DAG node to the given object's
  /// particle DAG node with an L->T edge.
  ///
  /// \param source - the DAGInfo object containing the particle node in
  ///                 question
  void LtoT(DAGInfo *target) {
    assert(has_normal());
    assert(target->has_parts());
    link_nodes(this, normal_, target, target->parts_, Operation::LtoT);
  }

  /// Create an S->T link in the DAG
  ///
  /// This will connect the particle DAG node of the given object to this
  /// object's particle DAG node with an S->T edge.
  ///
  /// \param source - the DAGInfo object containing the particle node in
  ///                 question
  void StoT(DAGInfo *source) {
    assert(source->has_parts());
    assert(has_parts());
    link_nodes(source, source->parts_, this, parts_, Operation::StoT);
  }

  /// Create an M->I link in the DAG
  ///
  /// This will connect the normal DAG node of the given object to this
  /// object's intermediate DAG node with an M->I edge.
  ///
  /// \param source - the DAGInfo object containing the normal node in
  ///                 question
  void MtoI(DAGInfo *source) {
    assert(source->has_normal());
    assert(has_interm());
    link_nodes(source, source->normal_, this, interm_, Operation::MtoI);
  }

  /// Create an I->I link in the DAG
  ///
  /// This will connect the intermediate DAG node of the given object to this
  /// object's intermediate DAG node with an I->I edge.
  ///
  /// \param source - the DAGInfo object containing the intermediate node in
  ///                 question
  void ItoI(DAGInfo *source) {
    assert(has_interm());
    assert(source->has_interm());
    link_nodes(source, source->interm_, this, interm_, Operation::ItoI);
  }

  /// Create an I->L link in the DAG
  ///
  /// This will connect the intermediate DAG node of the given object to this
  /// object's normal DAG node with an I->L edge.
  ///
  /// \param source - the DAGInfo object containing the intermediate node in
  ///                 question
  void ItoL(DAGInfo *source) {
    assert(source->has_interm());
    assert(has_normal());
    link_nodes(source, source->interm_, this, normal_, Operation::ItoL);
  }

  /// Collect the DAG nodes from this object
  ///
  /// During the realization of the DAG as LCO objects, the nodes are collected
  /// from each DAGInfo object. This routine serves that collection.
  ///
  /// \param terminals - the vector containing the source or target nodes
  /// \param internals - the vector containing the rest of the nodes
  void collect_DAG_nodes(std::vector<DAGNode *> &terminals,
                         std::vector<DAGNode *> &internals) {
    if (normal_) {
      internals.push_back(normal_);
    }
    if (interm_) {
      internals.push_back(interm_);
    }
    if (parts_) {
      terminals.push_back(parts_);
    }
  }

  /// Utility routine to connect nodes of the DAG
  ///
  /// This connects the given node with the specified operation.
  ///
  /// \param src_info - the DAGInfo object containing the source node
  /// \param source - the DAGNode that is the source of the edge being added
  /// \param dest_info - the DAGInfo object containing the destination node
  /// \param dest - the DAGNode that is the destination of the edge being added
  /// \param op - the operation to perform along the edge
  static void link_nodes(DAGInfo *src_info, DAGNode *source,
                         DAGInfo *dest_info, DAGNode *dest, Operation op) {
    dest_info->lock();
    dest->add_in_edge(source, op);
    dest_info->unlock();

    src_info->lock();
    source->add_out_edge(dest, op);
    src_info->unlock();
  }

 private:
  void lock() {
    hpx_lco_sema_p(lock_);
  }

  void unlock() {
    hpx_lco_sema_v(lock_, HPX_NULL);
  }

  Index idx_;
  hpx_addr_t lock_;
  DAGNode *normal_;
  DAGNode *interm_;
  DAGNode *parts_;   // source or target
};


} // namespace dashmm


#endif // __DASHMM_DAG_INFO_H__
