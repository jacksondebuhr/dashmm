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
  BHMethod(double theta) : theta_{theta} { }

  int type() const override {return kMethodBH;}
  MethodSerialPtr serialize() const override;

  //BHMethod has no special requirements
  bool compatible_with(const ExpansionRef expand) const override {return true;}

  void generate(SourceNode &curr, const ExpansionRef expand) const override;
  void aggregate(SourceNode &curr, const ExpansionRef expand) const override;
  void inherit(TargetNode &curr, const ExpansionRef expand,
                       size_t which_child) const override { }
  void process(TargetNode &curr, std::vector<SourceNode *> &consider,
                       bool curr_is_leaf) const override;

  //BHMethod always calls for refinement
  bool refine_test(bool same_sources_and_targets, const TargetNode &curr,
                   const std::vector<SourceNode *> &consider) const override {
    return true;
  }

  //Note that this is pretty weird. We pass a reference to a reference.
  // This is because MAC may cause a local version of the expansion to be
  // created, and so we do not want to pass the ExpansionRef by value and
  // create that local copy twice.
  bool MAC(ExpansionRef &expand, double size, Point pos) const;
  Point nearest(Point scenter, Point tcenter, double tsize) const;

 private:
  double theta_;
}


} // namespace dashmm


#endif // __DASHMM_BH_METHOD_H__
