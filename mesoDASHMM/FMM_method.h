#ifndef __FMM_METHOD_H__
#define __FMM_METHOD_H__


#include <vector>

#include "expansion.h"
#include "method.h"
#include "node.h"
#include "point.h"


namespace dashmm {


class FMM_Method : public Method {
 public:
  void generate(SourceNode *curr, const Expansion *expand) const override;
  void aggregate(SourceNode *curr, const Expansion *expand) const override;
  void inherit(TargetNode *curr, const Expansion *proto,
               size_t which_child) const override;
  void process(TargetNode *curr, std::vector<SourceNode *> &consider,
               bool curr_is_leaf) const override;
  bool refine_test(bool same_sources_and_targets, const TargetNode *curr,
                   const std::vector<SourceNode *> &consider) const override;

  bool well_sep_test_asymmetric(Index smaller, Index larger) const;
  bool well_sep_test(Index source, Index target) const;
};


} //namespace dashmm


#endif
