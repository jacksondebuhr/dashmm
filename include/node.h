#ifndef __DASHMM_NODE_H__
#define __DASHMM_NODE_H__


#include "include/domaingeometry.h"
#include "include/expansion.h"
#include "include/index.h"
#include "include/particle.h"


namespace dashmm {


//This is annoying, but what can you do. We will process what actually makes
// in into the documentation anyway...
struct NodeData {
  DomainGeometry root_geo;
  Index idx;
  hpx_addr_t parent;
  hpx_addr_t child[8];

  //TODO: Is this just a future that stores the ObjectBase?
  // then once it is ready, the future is set, and those reading it
  // can go ahead and get the information? I think so...
  hpx_addr_t expansion;
  hpx_addr_t method;

  //For Source Nodes, these refer to the source points
  // For target nodes, these refer to the target points
  hpx_addr_t parts;
  int n_parts;
  int n_parts_total;
};



class SourceNode {
 public:
  SourceNode(DomainGeometry g, int ix, int iy, int iz, int level,
             hpx_addr_t met, SourceNode *parent);
  explicit SourceNode(hpx_addr_t data = HPX_NULL)
      : local_{nullptr}, data_{data} { }
  ~SourceNode();

  void destroy();

  hpx_addr_t data() const {return data_;}

  //all the action is here
  // this returns an LCO to signal that partition is complete (this does not
  // mean that generation and aggregation are complete)
  hpx_addr_t partition(hpx_addr_t parts, int n_parts, int limit,
                       hpx_addr_t expand);

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

  mutable NodeData *local_;

  hpx_addr_t data_;
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
