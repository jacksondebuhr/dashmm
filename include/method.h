#ifndef __DASHMM_METHOD_H__
#define __DASHMM_METHOD_H__


#include <vector>


namespace dashmm {


//Probably these actions will target the node in question in most cases
//
//Or do we make these functions, and then the method actions will call these
// particular functions? What is the easiest thing when users start registering
// their own methods? I think we do not want the users to have to write actions.
// This way they are only having to do an extremely minimal amount of HPX-5
// if needed. Really we will likely go ahead and have a translator that will
// take some code markups and add what is needed to then to hide the HPX-5
// for defining the functions.
//
//Further for those that take node references, what do we do about that?
// have the functions written using the class? Probably...
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
  //Generall speaking, built-in methods will not be constructed this way by
  // the user. This form should only be used to create the first instance of
  // a method. For the built-in methods, this is handled by dashmm. For user-
  // defined methods, the user will need to perform this at least once.
  //
  //Though, this will also occur when we copy a method remotely. So we need
  // some mechanism to do that easily...
  Method(int type, std::vector<double> parms);
  // Otherwise, we generally just do copy construction from a prototype.

  bool valid() const;

  hpx_addr_t compatible_with(hpx_addr_t sync, const Expansion *expand) const;
  hpx_addr_t generate(hpx_addr_t sync, SourceNode *curr,
                      const Expansion *expand) const;
  hpx_addr_t aggregate(hpx_addr_t sync, SourceNode *curr,
                       const Expansion *expand) const;
  hpx_addr_t inherit(hpx_addr_t sync, TargetNode *curr, const Expansion *expand,
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


//TODO: consider moving this into an advanced interface file, or something
// like that. Perhaps something more like builtins?
int register_method(size_t params, hpx_action_t compat,
                    hpx_action_t gen, hpx_action_t agg, hpx_action_t inherit,
                    hpx_action_t proc, hpx_action_t reftest);


} // namespace dashmm


#endif // __DASHMM_METHOD_H__
