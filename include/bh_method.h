#ifndef __DASHMM_BH_METHOD_H__
#define __DASHMM_BH_METHOD_H__


#include "include/builtin_ids.h"
#include "include/expansionref.h"
#include "include/method.h"
#include "include/node.h"


namespace dashmm {


class BHMethod {
 public:
  //
  BHMethod(double theta);

  int type() const {return kMethodBH;}
  MethodSerialPtr serialize() const;

  //BHMethod has no special requirements
  bool compatible_with(const ExpansionRef expand) const {return true;}

  void generate(SourceNode &curr, const ExpansionRef expand) const;
  void aggregate(SourceNode &curr, const ExpansionRef expand) const;
  void inherit(TargetNode &curr, const ExpansionRef expand,
                       size_t which_child) const;
  void process(TargetNode &curr, std::vector<SourceNode *> &consider,
                       bool curr_is_leaf) const;

  //BHMethod always calls for refinement
  bool refine_test(bool same_sources_and_targets,
                           const TargetNode &curr,
                           const std::vector<SourceNode *> &consider) const {
    return true;
  }

 private:
  //
}


} // namespace dashmm


#endif // __DASHMM_BH_METHOD_H__
