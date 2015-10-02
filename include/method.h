#ifndef __DASHMM_METHOD_H__
#define __DASHMM_METHOD_H__


#include <vector>


namespace dashmm {


struct MethodDesc {
  size_t param_count;
  hpx_action_t compatible_with_function;
  hpx_action_t generate_function;
  hpx_action_t aggregate_function;
  hpx_action_t inherit_function;
  hpx_action_t process_function;
  hpx_action_t refine_test_function;
};


class Method {
 public:
  //Types that the functions implementing these need to take
  typedef bool (*compatible_with_handler_t)(const Expansion &);
  typedef void (*generate_handler_t)(SourceNode &, const Expansion &);
  typedef void (*aggregate_handler_t)(SourceNode &, const Expansion &);
  typedef void (*inherit_handler_t)(TargetNode &, const Expansion &, int which_child);

  //Generall speaking, built-in methods will not be constructed this way by
  // the user. This form should only be used to create the first instance of
  // a method. For the built-in methods, this is handled by dashmm. For user-
  // defined methods, the user will need to perform this at least once.
  //
  //Though, this will also occur when we copy a method remotely. So we need
  // some mechanism to do that easily...
  Method(int type, std::vector<double> parms);
  // Otherwise, we generally just do copy construction from a prototype.

  bool valid() const {return type != kFirstMethodType;}
  int type() const {return type_;}

  bool compatible_with(const Expansion &expand) const;
  hpx_addr_t generate(hpx_addr_t sync, SourceNode &curr,
                      const Expansion &expand) const;
  hpx_addr_t aggregate(hpx_addr_t sync, SourceNode &curr,
                       const Expansion &expand) const;
  hpx_addr_t inherit(hpx_addr_t sync, TargetNode &curr, const Expansion &expand,
                     size_t which_child) const;

  hpx_addr_t process(hpx_addr_t sync, TargetNode *curr,
                     std::vector<SourceNode *> &consider,
                     bool curr_is_leaf) const;
  hpx_addr_t refine_test(hpx_addr_t sync, bool same_sources_and_targets,
                           const TargetNode *curr,
                           const std::vector<SourceNode *> &consider) const;

 private:
  int type_;
  const MethodDesc &table_;
  std::vector<double> params_;
};


int register_method(MethodDesc desc);


} // namespace dashmm


#endif // __DASHMM_METHOD_H__
