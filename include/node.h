#ifndef __DASHMM_NODE_H__
#define __DASHMM_NODE_H__


/// \file include/node.h
/// \brief Source and Target trees


#include "include/domaingeometry.h"
#include "include/expansion.h"
#include "include/index.h"
#include "include/particle.h"


namespace dashmm {


/// The data for source nodes.
///
/// The data stored for SourceNode objects. This will be saved in a block in
/// the GAS.
struct SourceNodeData {
  /// The geometry of the root for the tree of which this node is a part.
  DomainGeometry root_geo;
  /// The index giving which subdivision of the root this node is.
  Index idx;
  /// The global address of the parent of this node.
  hpx_addr_t parent;
  /// The global addresses of the children of this node.
  hpx_addr_t child[8];
  /// The global address of the expansion object for this node.
  hpx_addr_t expansion;
  /// The global address of the method object for this node.
  hpx_addr_t method;
  /// A reference to the sources for this node. If this is an internal node,
  /// this will be an invalid SourceRef.
  SourceRef sources;
};


/// A node of the source tree.
///
/// This object is a reference object, storing only the HPX address of the
/// data for this node. As needed, the methods in this object will pin the
/// global data to produce a locally accessible version of the data to return
/// results in queries. As it is a reference, it should be passed by value
/// wherever it is needed.
class SourceNode {
 public:
  /// Construct from basic node data.
  SourceNode(DomainGeometry g, int ix, int iy, int iz, int level,
             hpx_addr_t met, SourceNode *parent);

  /// Construct from a global address.
  explicit SourceNode(hpx_addr_t data = HPX_NULL)
      : local_{nullptr}, data_{data} { }

  /// Copy constructor
  ///
  /// This is needed to assure that pin()/unpin() are properly reference counted
  /// when new instances of SourceNode refer to the same global data. If
  /// local_ was accidentally copied from a pinned node, the unpin() on either
  /// would invalidate the local pointer for the other.
  SourceNode(const SourceNode &other)
      : local_{nullptr}, data_{other.data()} { }

  /// Destructor
  ///
  /// This is required because the object may have pinned a global address,
  /// and so this address will have to be unpinned.
  ~SourceNode();

  /// Destroy the referred to data
  ///
  /// To delete the node, one must call destroy(); the destruction of the
  /// SourceNode object only destroys the reference to the global data, not
  /// the global data itself.
  void destroy();

  /// Return the address of the referred to data.
  hpx_addr_t data() const {return data_;}

  //all the action is here
  // this returns an LCO to signal that partition is complete (this does not
  // mean that generation and aggregation are complete)

  /// Partition the node.
  ///
  /// In addition to partitioning the domain based on the input particles, and
  /// input partition limit, this also schedules the work based on the method
  /// and the given expansion. In particular, the generate() and aggregate()
  /// methods will be called for this node based on the leaf/internal status
  /// of this node after partitioning.
  ///
  /// Partitioning will finish before the work of the method is finished, and
  /// for internal nodes, the difference in completion time could be large
  /// because the child node may have to wait itself for moments from its
  /// children before it can aggregate its own moments.
  ///
  /// This will return the address of an LCO that represents completion of the
  /// partitioning of this node.
  ///
  /// TODO: Is this correct?
  hpx_addr_t partition(hpx_addr_t parts, int n_parts, int limit,
                       hpx_addr_t expand);

  /// Is the referred to node a leaf?
  bool is_leaf() const;

  /// Does this refer to data?
  bool is_valid() const {return data_ != HPX_NULL;}

  DomainGeometry root_geo() const;
  Index index() const;
  int x_index() const;
  int y_index() const;
  int z_index() const;
  int level() const;
  SourceNode child(size_t i) const;
  SourceNode parent() const;
  ExpansionRef expansion() const;
  SourceRef parts() const;
  int n_parts() const;
  int n_parts_total() const;
  Point low() const;
  Point high() const;
  Point center() const;
  double size() const;

  //mutators
  void set_expansion(std::unique_ptr<Expansion> expand);

 private:
  void pin() const;
  void unpin() const;

  mutable SourceNodeData *local_;

  hpx_addr_t data_;
};


/// The data for target nodes.
///
/// The data stored for TargetNode objects. This will be saved in a block in
/// the GAS.
struct TargetNodeData {
  /// The geometry of the root for the tree of which this node is a part.
  DomainGeometry root_geo;
  /// The index giving which subdivision of the root this node is.
  Index idx;
  /// The global address of the parent of this node.
  hpx_addr_t parent;
  /// The global addresses of the children of this node.
  hpx_addr_t child[8];
  /// The global address of the expansion object for this node.
  hpx_addr_t expansion;
  /// The global address of the method object for this node.
  hpx_addr_t method;
  /// A reference to the targets for this node.
  TargetRef targets;
};


class TargetNode {
 public:
  TargetNode(DomainGeometry g, int ix, int iy, int iz, int level,
             hpx_addr_t method, TargetNode *parent);
  explicit TargetNode(hpx_addr_t data)
      : local_{nullptr}, data_{data} { }
  ~TargetNode();

  void destroy();

  hpx_addr_t data() const {return data_;}

  //all the action is here
  void partition(hpx_addr_t parts, int n_parts, int limit,
                 hpx_addr_t expand, int which_child,
                 bool same_sources_and_targets,
                 std::vector<SourceNode> consider);

  //queries
  bool is_leaf() const;
  bool is_valid() const {return data_ != HPX_NULL;}

  //accessors
  DomainGeometry root_geo() const;
  Index index() const;
  int x_index() const;
  int y_index() const;
  int z_index() const;
  int level() const;
  TargetNode child(size_t i) const;
  TargetNode parent() const;
  ExpansionRef expansion() const;
  TargetRef parts() const;
  int n_parts() const;
  int n_parts_total() const;
  Point low() const;
  Point high() const;
  Point center() const;
  double size() const;

  //mutators
  void set_expansion(std::unique_ptr<Expansion> expand);

 private:
  void pin() const;
  void unpin() const;

  mutable NodeData *local_;

  hpx_addr_t data_;
};


} // namespace dashmm


#endif // __DASHMM_NODE_H__
