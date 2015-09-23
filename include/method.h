#ifndef __DASHMM_METHOD_H__
#define __DASHMM_METHOD_H__


namespace dashmm {


class Method{
 public:
  virtual ~Method() { };

  // compatible with -
  //   This will check if the given expansion is compatible with the method.
  //   Typically this means a check to see if the required operators are
  //   given nontrivial implementations. This is signaled via the provides_L
  //   method of the expansion.
  virtual bool compatible_with(Expansion *expand) const = 0;

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


} // namespace dashmm


#endif // __DASHMM_METHOD_H__
