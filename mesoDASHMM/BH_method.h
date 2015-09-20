#ifndef __BH_METHOD_H__
#define __BH_METHOD_H__


#include "method.h"
#include "node.h"

#include <vector>


namespace dashmm {


class BH_Method : public Method {
 public:
  BH_Method(double crit);
  ~BH_Method();

  void generate(SourceNode *curr, const Expansion *expand) const override;
  void aggregate(SourceNode *curr, const Expansion *expand) const override;
  void inherit(TargetNode *curr, const Expansion *proto,
               size_t which_child) const override { }
  void process(TargetNode *curr, std::vector<SourceNode *> &consider,
               bool curr_is_leaf) const override;
  bool refine_test(bool same_sources_and_targets, const TargetNode *curr,
                   const std::vector<SourceNode *> &consider) const override {
    return true;
  };

  bool MAC(SourceNode *source, Point pos) const;

  Point nearest(const SourceNode *source, const TargetNode *target) const;

  void finalize(TargetNode *curr, std::vector<SourceNode *> &consider) const;

  double crit() const {return crit_;}

 private:
  double crit_;
};


} //namespace dashmm


#endif // __BH_METHOD_H__
