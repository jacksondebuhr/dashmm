#ifndef __DASHMM_NODE_H__
#define __DASHMM_NODE_H__


/// \file include/node.h
/// \brief Source and Target trees


#include "include/domaingeometry.h"
#include "include/expansion.h"
#include "include/expansionref.h"
#include "include/index.h"
#include "include/particle.h"
#include "include/types.h"


namespace dashmm {


// Forward declaration of internal types
struct SourceNodeData;
struct TargetNodeData;


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
  SourceNode(DomainGeometry g, Index idx,
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
  /// \param parts - the source points
  /// \param n_parts - the number of sources
  /// \param limit - the refinement limit for the partitioning
  /// \param type - the type of the expansion
  /// \param expand - global address to the prototype expansion data
  /// \param n_digits - accuracy of the expansion
  hpx_addr_t partition(Source *parts, int n_parts, int limit,
                       int type, hpx_addr_t expand, int n_digits);

  /// Is the referred to node a leaf?
  bool is_leaf() const;

  /// Does this refer to data?
  bool is_valid() const {return data_ != HPX_NULL;}

  /// Return the geometry of the root of the tree
  DomainGeometry root_geo() const;

  /// Return the index of this node
  Index index() const;

  /// Return the x index of this node
  int x_index() const;

  /// Return the y index of this node
  int y_index() const;

  /// Return the z index of this node
  int z_index() const;

  /// Return the level of this node
  int level() const;

  /// Return reference to a child of this node
  SourceNode child(size_t i) const;

  /// Return reference to the parent of this node
  SourceNode parent() const;

  /// Return this node's expansion
  ExpansionRef expansion() const;

  /// Return the source points for this node
  ///
  /// For a leaf node, this will return a reference to the source poitns in
  /// that node. For an internal node, this will return an empty reference
  SourceRef parts() const;

  /// Return the number of source points for this node
  ///
  /// For a lead node, this will return the number of source points in that
  /// leaf. For an internal node, this will return 0.
  int n_parts() const;

  /// Return the low corner of this node
  Point low() const;

  /// Return the high corner of this node
  Point high() const;

  /// Return the center of this node
  Point center() const;

  /// Return the edge size of this node
  double size() const;

  /// Set the expansion of this node.
  ///
  /// This sets the expansion of this node to be a reference to a globalized
  /// version of the probided expansion.
  ///
  /// \param expand - the expansion to set as the expansion of this node
  void set_expansion(std::unique_ptr<Expansion> expand);

 private:
  void pin() const;
  void unpin() const;

  mutable SourceNodeData *local_;

  hpx_addr_t data_;
};


class TargetNode {
 public:
  /// Construct from basic node data.
  ///
  /// This will create a node in the global address space with the specified
  /// data. Note that the expansion for the node will not be created; the
  /// expansion will have to be created in either inherit() or process().
  TargetNode(DomainGeometry g, Index idx, hpx_addr_t method,
             TargetNode *parent);

  /// Copy constructor
  ///
  /// This is needed so that there is correct reference counting for the
  /// pinned data locally.
  TargetNode(const TargetNode &other) : local_{nullptr}, data_{other.data()} { }

  /// Construct from a global address
  explicit TargetNode(hpx_addr_t data)
      : local_{nullptr}, data_{data} { }

  /// Destructor
  ///
  /// This is required as the object may have pinned a global address, so we
  /// need to allow for the opportunity to unpin the data.
  ~TargetNode();

  /// Destroy the referred data
  ///
  /// To delete the node, one must call destroy(); the destruction of the
  /// TargetNode object only destroys the reference to the global data, not
  /// the global data itself.
  void destroy();

  /// Return the address of the referred data.
  hpx_addr_t data() const {return data_;}

  /// Partition the node.
  ///
  /// In addition to partitioning the domain based on the input target
  /// locations, this also schedules the work based on the method for the
  /// various expansions, and the target points. In particular, inherit()
  /// and process() from the method associated with the node will be called.
  ///
  /// This call is synchronous. All contributions will have been made to all
  /// targets when this returns to the calling function.
  ///
  /// Note that for internal reasons, the target points are passed to this
  /// routine as a global address, rather than as the data itself.
  ///
  /// \param parts - global address of the target data
  /// \param n_parts - the number of target locations
  /// \param limit - the refinement limit for the partition
  /// \param expand - the expansion to use as prototype for the nodes created
  ///                 during the partition
  /// \param which_child - which child of this node's parent was this node
  /// \param same_sources_and_targets - are we in the case where the sources
  ///                 and targets are the same?
  /// \param consider - the list of source nodes in the consider list
  void partition(hpx_addr_t parts, int n_parts, int limit,
                 ExpansionRef expand, int which_child,
                 bool same_sources_and_targets,
                 std::vector<SourceNode> consider);

  /// Is the node a leaf?
  bool is_leaf() const;

  /// Is the reference referring to data?
  bool is_valid() const {return data_ != HPX_NULL;}

  /// Return the root geometry of this node's tree
  DomainGeometry root_geo() const;

  /// Return the index of this node
  Index index() const;

  /// Return the x index of this node
  int x_index() const;

  /// Return the y index of this node
  int y_index() const;

  /// Return the z index of this node
  int z_index() const;

  /// Return the level of this node
  int level() const;

  /// Return a reference to a child of this node
  TargetNode child(size_t i) const;

  /// Return a reference to the parent of this node
  TargetNode parent() const;

  /// Return a reference to this node's expansion
  ExpansionRef expansion() const;

  /// Return a reference to this node's target points.
  ///
  /// If this is an internal node, the returned reference will be invalid.
  /// If this is a leaf node, the returned reference will be valid.
  TargetRef parts() const;

  /// Return the number of targets in this node
  ///
  /// If this is a leaf node, this will return the number of targets contained
  /// in this node. Otherwise, this will return 0.
  int n_parts() const;

  /// Return the low corner of the node
  Point low() const;

  /// Return the high corner of the node
  Point high() const;

  /// Return the center of the node
  Point center() const;

  /// Returnd the side length of this node.
  double size() const;

  /// Set the expansion of this node.
  ///
  /// This will set this node's expansion to be a reference to the globalized
  /// form of the given expansion.
  void set_expansion(std::unique_ptr<Expansion> expand);

  /// Collect results into user's array
  void collect_results(hpx_addr_t user_array, size_t phi_offset);

 private:
  void pin() const;
  void unpin() const;

  mutable TargetNodeData *local_;

  hpx_addr_t data_;
};


} // namespace dashmm


#endif // __DASHMM_NODE_H__
