#ifndef __DASHMM_NODE_H__
#define __DASHMM_NODE_H__


#include "include/domaingeometry.h"
#include "include/expansion.h"
#include "include/index.h"
#include "include/particle.h"


namespace dashmm {


class SourceNode {
 public:
  SourceNode(DomainGeometry g, int ix, int iy, int iz, int level,
             hpx_addr_t met, SourceNode *parent);
  explicit SourceNode(hpx_addr_t data) : data_{data} { }
  ~SourceNode();

  hpx_addr_t data() const {return data_;}

  //all the action is here
  // this returns an LCO to signal that partition is complete (this does not
  // mean that generation and aggregation are complete)
  hpx_addr_t partition(hpx_addr_t first, hpx_addr_t last, int limit,
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
  SourceNode *child(size_t i) const;
  SourceNode *parent() const;
  hpx_addt_t expansion() const;
  hpx_addr_t first() const;
  hpx_addr_t last() const;
  Point low() const;
  Point high() const;
  Point center() const;
  double size() const;

  //mutators
  void set_expansion(hpx_addr_t expand);


 private:
  hpx_addr_t data_;
};


class TargetNode {
 public:
  TargetNode(DomainGeometry g, int ix, int iy, int iz, int level,
             hpx_addr_t method, TargetNode *parent);
  explicit TargetNode(hpx_addr_t data) : data_{data} { }
  ~TargetNode();

  hpx_addr_t data() const {return data_;}

  //all the action is here
  void partition(hpx_addr_t first, hpx_addr_t last, int limit,
                 hpx_addr_t expand, int which_child,
                 std::vector<SourceNode *> consider);

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
  TargetNode *child(size_t i) const;
  TargetNode *parent() const;
  hpx_addr_t expansion() const;
  hpx_addr_t first() const;
  hpx_addr_t last() const;
  Point low() const;
  Point high() const;
  Point center() const;
  double size() const;

  //mutators
  void set_expansion(hpx_addr_t expand);

 private:
  hpx_addr_t data_;
};


} // namespace dashmm


#endif // __DASHMM_NODE_H__
