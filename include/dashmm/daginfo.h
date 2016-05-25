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


#ifndef __DASHMM_DAG_INFO_H__
#define __DASHMM_DAG_INFO_H__


/// \file include/dashmm/daginfo.h
/// \brief Interface for intermediate representation of DAG


#include <vector>

#include <hpx/hpx.h>


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
};


struct DAGNode;

/// Edge in the explicit representation of the DAG
///
/// This is simply the target of the edge and the operation to peform on the
/// edge. The source of the edge is implicit in that the DAGNode which has
/// this edge in its list will be the source of the edge.
struct DAGEdge {
  const DAGNode *target;
  Operation op;

  DAGEdge() : target{nullptr}, op{Operation::Nop} { }
  DAGEdge(const DAGNode *end, Operation inop) : target{end}, op{inop} { }
};


/// Node in the explicit representation of the DAG
///
/// In addition to the set of edges that emanate from this node, the number
/// of inputs to this node is collected. This allows the eventual realization
/// of the DAG as ExpansionLCOs to have the correct number of inputs. Also,
/// when the distribution of the DAG is computed, the result will appear in
/// the locality entry of this object.
struct DAGNode {
  std::vector<DAGEdge> edges;    /// these are out edges
  int incoming;                  /// number of incoming edges
  int locality;                  /// the locality where this will be placed

  hpx_addr_t global_addx;        /// global address of object serving this node
  int other_member;              /// this is either n_digits for an expansion or
                                 /// n_targets for a target lco

  DAGNode()
      : edges{}, incoming{0}, locality{0}, global_addx{HPX_NULL},
        other_member{0} { }
  void add_input() {++incoming;}
  void add_edge(const DAGNode *end, Operation op) {
    edges.push_back(DAGEdge{end, op});
  }
};


/// DAG information relevant for a given tree node
///
/// Each node of the tree will have a DAGInfo member. This will store the
/// relevant information about the DAG for each tree node. This is the object
/// that Methods will use to construct the DAG. Each DAGInfo object will
/// represent up to three nodes of the DAG: one each for sources/targets,
/// the normal expansion, and an optional intermediate expansion.
///
/// A DAGInfo will create the normal expansion node when constructed. The
/// particle node or the intermediate node will be created on an as-needed
/// basis.
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
  DAGInfo() {
    normal_ = new DAGNode{};
    assert(normal_ != nullptr);
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
    delete normal_;
    normal_ = nullptr;
    if (interm_ != nullptr) {
      delete interm_;
      interm_ = nullptr;
    }
    if (parts_ != nullptr) {
      delete parts_;
      parts_ = nullptr;
    }
  }

  /// Add an intermediate node
  ///
  /// This will add an intermediate DAG node for the tree node associated
  /// with this object.
  void add_interm() {
    assert(interm_ == nullptr);
    interm_ = new DAGNode{};
    assert(interm_ != nullptr);
  }

  /// Add a particle node
  ///
  /// This will add a particle DAG node for the tree node associated with this
  /// object. This will represent either sources in the source tree, or
  /// targets in the target tree.
  void add_parts(int loc) {
    assert(parts_ == nullptr);
    parts_ = new DAGNode{};
    parts_->locality = loc;
    assert(parts_ != nullptr);
  }

  // TODO: be aware that these might need to change
  DAGInfo(const DAGInfo &other) = delete;
  DAGInfo &operator=(const DAGInfo &other) = delete;
  DAGInfo(const DAGInfo &&other) = delete;
  DAGInfo &operator=(const DAGInfo &&other) = delete;

  bool has_interm() const {return interm_ != nullptr;}
  bool has_parts() const {return parts_ != nullptr;}

  const DAGNode *normal() const {return normal_;}
  const DAGNode *interm() const {return interm_;}
  const DAGNode *parts() const {return parts_;}

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

  /// Create an S->M link in the DAG
  ///
  /// This will connect the source DAG node of the given object to this object's
  /// normal DAG node with an S->M edge.
  ///
  /// \param source - the DAGInfo object containing the source node in question
  void StoM(DAGInfo *source) {
    assert(source->has_parts());
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
    link_nodes(source, source->parts_, this, normal_, Operation::StoL);
  }

  /// Create an M->M link in the DAG
  ///
  /// This will connect the normal DAG node of the given object to this object's
  /// normal DAG node with an M->M edge.
  ///
  /// \param source - the DAGInfo object containing the normal node in question
  void MtoM(DAGInfo *source) {
    link_nodes(source, source->normal_, this, normal_, Operation::MtoM);
  }

  /// Create an M->L link in the DAG
  ///
  /// This will connect the normal DAG node of the given object to this object's
  /// normal DAG node with an M->L edge.
  ///
  /// \param source - the DAGInfo object containing the normal node in question
  void MtoL(DAGInfo *source) {
    link_nodes(source, source->normal_, this, normal_, Operation::MtoL);
  }

  /// Create an L->L link in the DAG
  ///
  /// This will connect the normal DAG node of the given object to this object's
  /// normal DAG node with an L->L edge.
  ///
  /// \param source - the DAGInfo object containing the normal node in question
  void LtoL(DAGInfo *source) {
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

  /// Collect the DAG nodes from this object
  ///
  /// During the realization of the DAG as LCO objects, the nodes are collected
  /// from each DAGInfo object. This routine serves that collection.
  ///
  /// \param terminals - the vector containing the source or target nodes
  /// \param internals - the vector containing the rest of the nodes
  void collect_DAG_nodes(std::vector<DAGNode *> &terminals,
                         std::vector<DAGNode *> &internals) {
    internals.push_back(normal_);
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
    dest->add_input();
    dest_info->unlock();

    src_info->lock();
    source->add_edge(dest, op);
    src_info->unlock();
  }

 private:
   void lock() {
     hpx_lco_sema_p(lock_);
   }

   void unlock() {
     // TODO: can I get away with nonsync and ignore sync on this?
     hpx_lco_sema_v_sync(lock_);
   }

  hpx_addr_t lock_;
  // We do these as heap allocations so that only those that are needed
  // are created.
  DAGNode *normal_;
  DAGNode *interm_;
  DAGNode *parts_;   // source or target
};


} // namespace dashmm


#endif // __DASHMM_DAG_INFO_H__
