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
  Nop;
  StoM;
  StoL;
  MtoM;
  MtoL;
  LtoL;
  MtoT;
  LtoT;
  StoT;
};


struct DAGNode;

/// Edge in the explicit representation of the DAG
///
/// This is simply the target of the edge and the operation to peform on the
/// edge. The source of the edge is implicit in that the DAGNode which has
/// this edge in its list will be the source of the edge.
struct DAGEdge {
  // TODO: Should this be const?
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
  std::vector<DAGEdge> edges;    // these are out edges
  int incoming;                  // number of incoming edges
  int locality;                  // the locality where this will be placed

  // This is less than ideal, but either object that will be stored only has
  // two members, one an hpx_addr_t and one an int. So this will be the solution
  // for now.
  hpx_addr_t global_addx;        // global address of object serving this node
  int other_member;              // this is either n_digits for an expansion or
                                 // n_targets for a target lco

  DAGNode()
      : edges{}, incoming{0}, locality{0}, global_addx{HPX_NULL},
        other_member{0} { }
  void add_input() {++incoming;}
  void add_edge(const DAGNode *end, Operation op) {
    edges.push_back(DAGEdge{end, op})
  }
};


class DAGInfo {
 public:
  // Construct as needed
  // NOTE: This can only be called from an HPX thread.
  DAGInfo() {
    normal_ = new DAGNode{};
    assert(normal_ != nullptr);
    expon_ = nullptr;
    parts_ = nullptr;
    lock_ = hpx_lco_sema_new(1);
    assert(lock_ != HPX_NULL);
  }


  ~DAGInfo() {
    hpx_lco_delete_sync(lock_);
    delete normal_;
    normal_ = nullptr;
    if (expon_ != nullptr) {
      delete expon_;
      expon_ = nullptr;
    }
    if (parts_ != nullptr) {
      delete parts_;
      parts_ = nullptr;
    }
  }

  void add_expon() {
    assert(expon_ == nullptr);
    expon_ = new DAGNode{};
    assert(expon_ != nullptr);
  }

  void add_parts() {
    assert(parts_ == nullptr);
    parts_ = new DAGNode{};
    assert(parts_ != nullptr);
  }

  // TODO: be aware that these might need to change
  // Copy construction and assignment are removed
  DAGInfo(const DAGInfo &other) = delete;
  DAGInfo &operator=(const DAGInfo &other) = delete;
  // Move construction and assignment are removed
  DAGInfo(const DAGInfo &&other) = delete;
  DAGInfo &operator=(const DAGInfo &other) = delete;

  bool has_expon() const {return expon_ != nullptr;}
  bool has_parts() const {return parts_ != nullptr;}

  const DAGNode *normal() {return normal_;}
  const DAGNode *expon() {return expon_;}
  const DAGNode *parts() {return parts_;}

  void StoM(DAGInfo *source) {
    assert(source->has_parts());
    link_nodes(source->parts_, normal_, Operation::StoM);
  }

  void StoL(DAGInfo *source) {
    assert(source->has_parts());
    link_nodes(source->parts_, normal_, Operation::StoL);
  }

  void MtoM(DAGInfo *source) {
    link_nodes(source->normal_, normal_, Operation::MtoM);
  }

  void MtoL(DAGInfo *source) {
    link_nodes(source->normal_, normal_, Operation::MtoL);
  }

  void LtoL(DAGInfo *source) {
    link_nodes(source->normal_, normal_, Operation::LtoL);
  }

  void MtoT(DAGInfo *target) {
    assert(target->has_parts());
    link_nodes(normal_, target->parts_, Operation::MtoT);
  }

  void LtoT(DAGInfo *target) {
    assert(target->has_parts());
    link_nodes(normal_, target->parts_, Operation::LtoT);
  }

  void StoT(DAGInfo *source) {
    assert(source->has_parts());
    assert(has_parts());
    link_nodes(source->parts_, parts_, Operation::StoT);
  }

  static void link_nodes(DAGNode *source, DAGNode *dest, Operation op) {
    dest->add_input();
    source->add_edge(dest, op);
  }

 private:
  hpx_addr_t lock_;
  // We do these as heap allocations so that only those that are needed
  // are created.
  DAGNode *normal_;
  DAGNode *expon_;
  DAGNode *parts_;   // source or target
}


} // namespace dashmm


#endif // __DASHMM_DAG_INFO_H__
