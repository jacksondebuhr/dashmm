#ifndef __METHOD_H__
#define __METHOD_H__


#include <cmath>
#include <cstdio>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

#include "expansion.h"
#include "particle.h"
#include "point.h"
#include "node.h"


namespace dashmm {


//TODO:
// In the HPX version, what is this? In principle the object is just a
// specification of the actions to take. But likely what we want is a copy of
// this specification at each locality. Then there is just some known mapping
// from the specific instance to the local copy. Perhaps we have an identifier
// for each "named" object in dashmm.
//
// Then dashmm::init is responsible for getting these created or something.
// and then these are shared around the system.
//
// Then asking for the method or expansion would share the needed data across
// the system, and register the name with dashmm so that local access is quick.
//
// So we could have a dashmm handle type (basically some kind of integer)
// and there are maps at each locality. The first is a map from that handle
// to the type: either GAS object or one that has a local copy.
// Then one of two maps are used to get the actual object reference. One map
// stores the hpx_addr_t of the global object. The other stores the local
// address of the shared local copy of a globally named thing.
//
// Note, that this seems a bit awkward. We shall need to look into what we do
// with these handles to see how this will go.

// I guess another option here is to go ahead and just copy the values around.
// Not sure what sort of overhead that implies. We will already be moving
// a fair bit of data about, so a dual lookup of the values might just be
// long enough that it is not worth it. One could always query dashmm for the
// official version of the thing, but we could just pass the action ids or
// whatever around.


/// Base class from which all methods are derived.
class Method{
 public:
  virtual ~Method() { };

  // Generate -
  //   This is the work that happens at the leaf nodes of the source tree
  //   during the ascent. This is probably just the generation of the moments
  //   from the sources in that leaf, but it could be something else.
  virtual void generate(SourceNode *curr, const Expansion *expand) const = 0;

  // Aggregate -
  //   This is the work that happens at the internal nodes of the source
  //   tree during the ascent. This is probably just the combination of the
  //   moments of the argument's children into the moment for the argument.
  virtual void aggregate(SourceNode *curr, const Expansion *expand) const = 0;

  // Inherit -
  //   This is any work that will modify the current node's expansion based
  //   on information from the parent. e.g., in FMM, this is the local to local
  //   from parent to child.
  virtual void inherit(TargetNode *curr, const Expansion *proto,
                       size_t which_child) const = 0;

  // Process -
  //   This takes in a vector of source nodes to consider, and takes action
  //   based on each. This might be to apply some operation based on the
  //   node, or to pass the node onto the children. The consider parameter
  //   will be modified by this function.
  virtual void process(TargetNode *curr, std::vector<SourceNode *> &consider,
                       bool curr_is_leaf) const = 0;

  // Refine Test -
  //
  virtual bool refine_test(bool same_sources_and_targets,
                           const TargetNode *curr,
                           const std::vector<SourceNode *> &consider) const = 0;
};


} //namespace dashmm


#endif // __METHOD_H__
