#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

#include "node.h"
#include "method.h"


namespace dashmm {


//*******************************************************************
// First the SourceNode implementation
//*******************************************************************


SourceNode::SourceNode(DomainGeometry g, int ix, int iy, int iz, int level,
                       const Method *m, SourceNode *parent)
    : root_geo_{g}, idx_{ix, iy, iz, level},
      parent_{parent}, child_{}, exp_{}, first_{}, last_{}, met_{m} { }


SourceNode::~SourceNode(){
  for (size_t i = 0; i < 8; ++i) {
    if (child_[i]) {
      delete child_[i];
      child_[i] = nullptr;
    }
  }
}


//The provided iterators do not represent a sorted sequence. They merely
// represent a sequence of particles in the volume of this node.
void SourceNode::partition(std::vector<Source>::iterator first,
                           std::vector<Source>::iterator last, int limit,
                           const Expansion *expand) {
  assert(last - first > 0);
  first_ = first;
  last_ = last;

  if ((last - first) <= limit) {
    //We need to generate moments here
    met_->generate(this, expand);
    return;
  }

  //first partition the particles among the children of this node
  std::vector<Source>::iterator splits[9]{};
  splits[0] = first;
  splits[8] = last;

  Point cen{center()};

  // first by z
  double z_center = cen.z();
  auto z_comp = [&z_center](Source a) {
    return a.z() < z_center;
  };
  splits[4] = std::partition(first, last, z_comp);

  // then by y
  double y_center = cen.y();
  auto y_comp = [&y_center](Source a) {
    return a.y() < y_center;
  };
  splits[2] = std::partition(first, splits[4], y_comp);
  splits[6] = std::partition(splits[4], last, y_comp);

  // then by x
  double x_center = cen.x();
  auto x_comp = [&x_center](Source a) {
    return a.x() < x_center;
  };
  splits[1] = std::partition(first, splits[2], x_comp);
  splits[3] = std::partition(splits[2], splits[4], x_comp);
  splits[5] = std::partition(splits[4], splits[6], x_comp);
  splits[7] = std::partition(splits[6], last, x_comp);

  //then create any needed children of this node
  for (size_t i = 0; i < 8; ++i) {
    //create a child if there is a non-empty sequence for that child
    if (splits[i] != splits[i + 1]) {
      Index cidx{idx_.child(i)};
      child_[i] = new SourceNode{root_geo_, cidx.x(), cidx.y(), cidx.z(),
                           cidx.level(), met_, this};
      child_[i]->partition(splits[i], splits[i + 1], limit, expand);
    } else {
      child_[i] = nullptr;
    }
  }

  assert(!is_leaf());
  met_->aggregate(this, expand);
}


bool SourceNode::is_leaf() const {
  for (size_t i = 0; i < 8; ++i) {
    if (child_[i]) {
      return false;
    }
  }
  return true;
}


//*******************************************************************
// Second the TargetNode implementation
//*******************************************************************


TargetNode::TargetNode(DomainGeometry g, int ix, int iy, int iz, int level,
                       const Method *m, TargetNode *parent)
    : root_geo_{g}, idx_{ix, iy, iz, level},
      parent_{parent}, child_{}, exp_{}, first_{}, last_{}, met_{m} { }


TargetNode::~TargetNode(){
  for (size_t i = 0; i < 8; ++i) {
    if (child_[i]) {
      delete child_[i];
      child_[i] = nullptr;
    }
  }
}


//The provided iterators do not represent a sorted sequence. They merely
// represent a sequence of particles in the volume of this node.
//TODO: In the hpx version, we will merge this with descend. Also, we will
// add a second criterion to the partition. If the neighbors in the source
// tree of this node are not refined, we will not refine here.
void TargetNode::partition(std::vector<Target>::iterator first,
                           std::vector<Target>::iterator last, int limit,
                           const Expansion *expand, int which_child,
                           std::vector<SourceNode *> consider) {
  assert(last - first > 0);
  first_ = first;
  last_ = last;

  bool refine = false;
  if ((last - first) > limit) {
    //TODO: fix that false
    refine = met_->refine_test(false, this, consider);
  }

  met_->inherit(this, expand, which_child);
  //If we are not going to refine, then curr is a leaf
  met_->process(this, consider, !refine);

  if (!refine) {
    return;
  }

  //first partition the particles among the children of this node
  std::vector<Target>::iterator splits[9]{};
  splits[0] = first;
  splits[8] = last;

  Point cen{center()};

  // first by z
  double z_center = cen.z();
  auto z_comp = [&z_center](Target a) {
    return a.z() < z_center;
  };
  splits[4] = std::partition(first, last, z_comp);

  // then by y
  double y_center = cen.y();
  auto y_comp = [&y_center](Target a) {
    return a.y() < y_center;
  };
  splits[2] = std::partition(first, splits[4], y_comp);
  splits[6] = std::partition(splits[4], last, y_comp);

  // then by x
  double x_center = cen.x();
  auto x_comp = [&x_center](Target a) {
    return a.x() < x_center;
  };
  splits[1] = std::partition(first, splits[2], x_comp);
  splits[3] = std::partition(splits[2], splits[4], x_comp);
  splits[5] = std::partition(splits[4], splits[6], x_comp);
  splits[7] = std::partition(splits[6], last, x_comp);

  //then create any needed children of this node
  for (size_t i = 0; i < 8; ++i) {
    //create a child if there is a non-empty sequence for that child
    if (splits[i] != splits[i + 1]) {
      Index cidx{idx_.child(i)};
      child_[i] = new TargetNode{root_geo_, cidx.x(), cidx.y(), cidx.z(),
                           cidx.level(), met_, this};
      child_[i]->partition(splits[i], splits[i + 1], limit, expand, i,
                           consider);
    } else {
      child_[i] = nullptr;
    }
  }
}


/*
void TargetNode::descend(const Expansion *proto, size_t which_child,
                   std::vector<SourceNode *> consider) {
  met_->inherit(this, proto, which_child);
  met_->process(this, consider);

  for (size_t i = 0; i < 8; ++i) {
    if (child_[i]) {
      child_[i]->descend(proto, i, consider);
    }
  }
}
*/


bool TargetNode::is_leaf() const {
  for (size_t i = 0; i < 8; ++i) {
    if (child_[i]) {
      return false;
    }
  }
  return true;
}


} //namespace dashmm
