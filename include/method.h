#ifndef __DASHMM_METHOD_H__
#define __DASHMM_METHOD_H__


/// \file include/method.h
/// \brief Abstract interface for Method objects


#include <vector>


#include "include/expansionref.h"
#include "include/node.h"
#include "include/types.h"


namespace dashmm {


/// The smallest allowed method type identifier for user-defined methods
extern constexpr int kFirstUserMethodType;

/// The largest allowed method type identifier for user-defined methods
extern constexpr int kLastUserMethodType;


/// The type for a serialized method object
///
/// Methods are all required to produce a serialized version of themselves
/// on request from the user. This serialized form will be of this type.
struct MethodSerial {
  int reserved;
  int type;
  size_t size;      //NOTE: This refers to the data, not the overall size here
  double data[];
};


/// The signature of a function that creates a method
///
/// To register a method with DASHMM, the user will need to provide a function
/// that can be used to create that sort of Method given data in the form of
/// a MethodSerial pointer.
typedef Method *(*method_creation_function_t)(size_t, MethodSerial *);


/// Abstract interface for Methods used in DASHMM
///
/// This interface specifies the requirements for methods that a user
/// may add to DASHMM. In general, the user will need to implement the method
/// and then somewhere in their program call register_method() providing the
/// creation functions that the user also implements. For details, see the
/// DASHMM Advanced User Guide.
class Method {
 public:
  virtual ~Method() { }

  /// Release the internal data for a method
  ///
  /// Method objects need to support the ability to export the data making up
  /// method in a contiguous chunk of memory. These serialized versions are
  /// used throughout the system. When release() is called, the method object
  /// is then 'empty' in a sense, and the method no longer owns any data.
  virtual MethodSerial *release() const = 0;

  /// The type identifier of the method
  virtual int type() const = 0;

  /// Determine if this method is compatible with a given expansion
  ///
  /// This will determine the compatability of the method with the given
  /// expansion. This is generally down to deciding if the method uses
  /// local or exponential operators, and checking that the expansion
  /// can provide those operators.
  ///
  /// \param expand - a reference to the expansion in question
  ///
  /// \returns - true in case of compatability; false otherwise
  virtual bool compatible_with(const ExpansionRef expand) const = 0;

  /// Generate the expansion at the leaf of the source tree
  ///
  /// This operation is invoked at the leaves of the source tree to generate
  /// the expansions for the sources. Typically, these will be multipole
  /// expansions.
  ///
  /// It is assumed that generate will either create the expansion with the
  /// needed data, or it will create an empty expansion and schedule the
  /// contribution to the expansion. Internally, DASHMM will call finalize on
  /// the expansion after generate() is called for a particular node.
  /// Further, generate will not return until the expansion for @p curr has
  /// been set. This does not require that all contributions have been made to
  /// that expansion.
  ///
  /// \param curr - the current node of the source tree (will be a leaf)
  /// \param expand - a reference to a prototype expansion that can be used
  ///                 to generate the expansion for the given node.
  virtual void generate(SourceNode &curr, const ExpansionRef expand) const = 0;

  /// Combine expansions from children of an internal source node
  ///
  /// This operation is invoked on internal nodes of the source tree to
  /// combine the expansions of the children on the given node into the
  /// expansion for the internal node.
  ///
  /// It is assumed that aggregate creates the expansion and schedules all
  /// contributions to that expansion.
  ///
  /// \param curr - the current node of the source tree (will be internal)
  /// \param expand - a reference to a prototype expansion that can be used
  ///                 to aggregate the expansions of the given node's children.
  virtual void aggregate(SourceNode &curr, const ExpansionRef expand) const = 0;

  /// Inherit an expansion from a target node's parent
  ///
  /// This operation is invoked on nodes in the target tree to inherit the
  /// effect of the expansion collected at the parent of the given node.
  ///
  /// \param curr - the current node of the target tree
  /// \param expand - a reference to a prototype expansion that might be used
  /// \which_child - which child @p curr is of its parent
  virtual void inherit(TargetNode &curr, const ExpansionRef expand,
                       size_t which_child) const = 0;

  /// Process the list of source nodes for a given target node
  ///
  /// The bulk of the work in the traversal of the target tree is in process.
  /// This operation is called on each target node and it will look at the
  /// set of source nodes in @p consider, and take the appropriate action
  /// given the usability of those nodes.
  ///
  /// \param curr - the target node in question
  /// \param consider - a vector of the source nodes under consideration
  /// \param curr_is_leaf - indicates if @p curr is a leaf node
  virtual void process(TargetNode &curr, std::vector<SourceNode> &consider,
                       bool curr_is_leaf) const = 0;

  /// Decide if the node in the target tree should be refined
  ///
  /// In some advanced methods, it can be desirable to halt refinement of the
  /// target tree before the specified refinement limit is reached. This
  /// operation performs that test.
  ///
  /// \param same_sources_and_targets - are the source and target points the
  ///                                   same
  /// \param curr - the target node being examined
  /// \param consider - the list of source nodes being considered at this node
  ///
  /// \returns - true if the refinement should proceed; false otherwise
  virtual bool refine_test(bool same_sources_and_targets,
                           const TargetNode &curr,
                           const std::vector<SourceNode> &consider) const = 0;
};


/// Register a user-defined method with DASHMM
///
/// This will connect the specified creation function @p creator with the
/// @p type for DASHMM. This will allow the system to interact with user-defined
/// methods. For details, see the DASHMM Advanced User Guide.
///
/// \param type - the type identifier of the method being registered.
/// \param creator - the action identifier of the creation function for the
///                  method being registered.
///
/// \returns - kSuccess on success; kDomainError otherwise, indicating that the
///            specified @p type is either out of range, or already in use.
ReturnCode register_method(int type, hpx_action_t creator);


//NOTE: Not intended for end-user use


/// Initialize the method registration table
void init_method_table();


/// Finalize (clean up) the method registration table
void fini_method_table()


/// Create a method object of the given type
///
/// This will create a new method object of the given type. The returned
/// method will be constructed from the create function provided by the user
/// during method registration.
///
/// \param type - the type of method
/// \param data - the serialized data around which the method will be created
///
/// \returns - the resulting method object
std::unique_ptr<Method> create_method(int type, MethodSerial *data);


} // namespace dashmm


#endif // __DASHMM_METHOD_H__
