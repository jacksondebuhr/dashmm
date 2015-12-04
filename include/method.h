#ifndef __DASHMM_METHOD_H__
#define __DASHMM_METHOD_H__


#include <vector>


#include "include/expansionref.h"
#include "include/node.h"
#include "include/types.h"


namespace dashmm {


extern constexpr int kFirstUserMethodType;
extern constexpr int kLastUserMethodType;


typedef Method *(*method_creation_function_t)(size_t, void *);


struct MethodSerial {
  int reserved;
  int type;
  size_t size;      //NOTE: This refers to the data, not the overall size here
  char data[];
};


class Method {
 public:
  virtual ~Method() { }

  virtual int type() const = 0;
  virtual MethodSerial *release() const = 0;

  virtual bool compatible_with(const ExpansionRef expand) const = 0;
  virtual void generate(SourceNode &curr, const ExpansionRef expand) const = 0;
  virtual void aggregate(SourceNode &curr, const ExpansionRef expand) const = 0;
  virtual void inherit(TargetNode &curr, const ExpansionRef expand,
                       size_t which_child) const = 0;
  virtual void process(TargetNode &curr, std::vector<SourceNode> &consider,
                       bool curr_is_leaf) const = 0;
  virtual bool refine_test(bool same_sources_and_targets,
                           const TargetNode &curr,
                           const std::vector<SourceNode> &consider) const = 0;
};


//returns true on success
ReturnCode register_method(int type, hpx_action_t creator);


//NOTE: Not intended for end-user use
std::unique_ptr<Method> create_method(int type, MethodSerial *data);


} // namespace dashmm


#endif // __DASHMM_METHOD_H__
