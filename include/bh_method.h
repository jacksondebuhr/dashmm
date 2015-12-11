#ifndef __DASHMM_BH_METHOD_H__
#define __DASHMM_BH_METHOD_H__


#include "include/ids.h"
#include "include/expansionref.h"
#include "include/method.h"
#include "include/node.h"


namespace dashmm {


class BHMethod : public Method {
 public:
  BHMethod(double theta) : local_{nullptr} {
    local_ = static_cast<MethodSerial *>(
            malloc(sizeof(MethodSerial) + sizeof(double)));
    assert(local_);
    local_->type = kMethodBH;
    local_->size = sizeof(double);
    local_->data[0] = theta;
  }

  MethodSerial *release() override {
    MethodSerial *retval{local_};
    local_ = nullptr;
    return retval;
  }

  size_t bytes() const {
    return (sizeof(MethodSerial) + sizeof(double));
  }

  int type() const override {return kMethodBH;}

  bool compatible_with(const Expansion *expand) const override {return true;}

  void generate(SourceNode &curr, const ExpansionRef expand) const override;
  void aggregate(SourceNode &curr, const ExpansionRef expand) const override;
  void inherit(TargetNode &curr, const ExpansionRef expand,
                       size_t which_child) const override { }
  void process(TargetNode &curr, std::vector<SourceNode> &consider,
                       bool curr_is_leaf) const override;

  //BHMethod always calls for refinement
  bool refine_test(bool same_sources_and_targets, const TargetNode &curr,
                   const std::vector<SourceNode> &consider) const override {
    return true;
  }

  //Note that this is pretty weird. We pass a reference to a reference.
  // This is because MAC may cause a local version of the expansion to be
  // created, and so we do not want to pass the ExpansionRef by value and
  // create that local copy twice.
  bool MAC(Point exp_point, double size, Point pos) const;
  Point nearest(Point scenter, Point tcenter, double tsize) const;

 private:
  MethodSerial *local_;
};


} // namespace dashmm


#endif // __DASHMM_BH_METHOD_H__
