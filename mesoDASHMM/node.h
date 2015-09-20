#ifndef __DASHMM_NODE_H__
#define __DASHMM_NODE_H__


#include <memory>
#include <vector>

#include "domaingeometry.h"
#include "expansion.h"
#include "index.h"
#include "particle.h"


namespace dashmm {


class Method;


class SourceNode{
 public:
  SourceNode(DomainGeometry g, int ix, int iy, int iz, int level,
       const Method *m, SourceNode *parent);
  ~SourceNode();

  void partition(std::vector<Source>::iterator first,
                 std::vector<Source>::iterator last, int limit,
                 const Expansion *expand);

  bool is_leaf() const;

  //accessors
  DomainGeometry root_geo() const {return root_geo_;}
  Index index() const {return idx_;}
  int x_index() const {return idx_.x();}
  int y_index() const {return idx_.y();}
  int z_index() const {return idx_.z();}
  int level() const {return idx_.level();}
  SourceNode *child(size_t i) const {return child_[i];}
  SourceNode *parent() const {return parent_;}
  Expansion *expansion() const {return exp_.get();}
  std::vector<Source>::iterator first() const {return first_;}
  std::vector<Source>::iterator last() const {return last_;}
  Point low() const {
    return root_geo_.low_from_index(idx_.x(), idx_.y(), idx_.z(), idx_.level());
  }
  Point high() const {
    return root_geo_.high_from_index(idx_.x(), idx_.y(), idx_.z(),
                                     idx_.level());
  }
  Point center() const {
    return root_geo_.center_from_index(idx_.x(), idx_.y(), idx_.z(),
                                       idx_.level());
  }
  double size() const {
    return root_geo_.size_from_level(idx_.level());
  }

  void set_expansion(std::unique_ptr<Expansion> expand) {
    exp_ = std::move(expand);
  }

 private:
  DomainGeometry root_geo_; //as is
  Index idx_;               //as is
  SourceNode *parent_;      //hpx_addr_t
  SourceNode *child_[8];    //hpx_addr_t
  std::unique_ptr<Expansion> exp_;  //TODO: What to do with this one?
  std::vector<Source>::iterator first_;   //hpx_addr_t or an index?
  std::vector<Source>::iterator last_;    //hpx_addr_t or an index?
  const Method *met_;               //TODO what here?
};


class TargetNode{
 public:
  TargetNode(DomainGeometry g, int ix, int iy, int iz, int level,
       const Method *m, TargetNode *parent);
  ~TargetNode();

  void partition(std::vector<Target>::iterator first,
                 std::vector<Target>::iterator last, int limit,
                 const Expansion *expand, int which_child,
                 std::vector<SourceNode *> consider);
  //void descend(const Expansion *proto, size_t which_child,
  //             std::vector<SourceNode *> consider);

  bool is_leaf() const;

  //accessors
  DomainGeometry root_geo() const {return root_geo_;}
  Index index() const {return idx_;}
  int x_index() const {return idx_.x();}
  int y_index() const {return idx_.y();}
  int z_index() const {return idx_.z();}
  int level() const {return idx_.level();}
  TargetNode *child(size_t i) const {return child_[i];}
  TargetNode *parent() const {return parent_;}
  Expansion *expansion() const {return exp_.get();}
  std::vector<Target>::iterator first() const {return first_;}
  std::vector<Target>::iterator last() const {return last_;}
  Point low() const {
    return root_geo_.low_from_index(idx_.x(), idx_.y(), idx_.z(), idx_.level());
  }
  Point high() const {
    return root_geo_.high_from_index(idx_.x(), idx_.y(), idx_.z(),
                                     idx_.level());
  }
  Point center() const {
    return root_geo_.center_from_index(idx_.x(), idx_.y(), idx_.z(),
                                       idx_.level());
  }
  double size() const {
    return root_geo_.size_from_level(idx_.level());
  }

  void set_expansion(std::unique_ptr<Expansion> expand) {
    exp_ = std::move(expand);
  }

 private:
  DomainGeometry root_geo_; //as is
  Index idx_;               //as is
  TargetNode *parent_;      //hpx_addr_t
  TargetNode *child_[8];    //hpx_addr_t
  std::unique_ptr<Expansion> exp_;  //TODO: What to do with this one?
  std::vector<Target>::iterator first_;   //hpx_addr_t or an index?
  std::vector<Target>::iterator last_;    //hpx_addr_t or an index?
  const Method *met_;               //TODO what here?
};


} //namespace dashmm


#endif // __DASHMM_NODE_H__
