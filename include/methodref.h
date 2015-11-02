#ifndef __DASHMM_METHOD_REF_H__
#define __DASHMM_METHOD_REF_H__


#include <vector>

#include <hpx/hpx.h>

#include "include/expansionref.h"
#include "include/method.h"
#include "include/node.h"


namespace dashmm {


class MethodRef {
 public:
  explicit MethodRef(hpx_addr_t addr) : data_{ref}, met_{nullptr} { }

  int type() const;
  bool compatible_with(const ExpansionRef expand) const;
  void generate(SourceNode &curr, const ExpansionRef expand) const;
  void aggregate(SourceNode &curr, const ExpansionRef expand) const;
  void inherit(TargetNode &curr, const ExpansionRef expand,
               size_t which_child) const;
  void process(TargetNode &curr, std::vector<SourceNode> &consider,
               bool curr_is_leaf) const;
  void bool refine_test(bool same_sources_and_targets,
                        const TargetNode &curr,
                        const std::vector<SourceNode> &consider) const;

 private:
  MethodSerial *pin();
  void unpin();
  void setup_local_method();

  hpx_addr_t data_;
  std::unique_ptr<Method> met_;
};


MethodRef globalize_method(Method *met, hpx_addr_t where);


} // namespace dashmm


#endif // __DASHMM_METHOD_REF_H__
