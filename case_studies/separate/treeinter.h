// This is the template interface thing to isolate tree stuff from DASHMM
//
// There are still some details to work out about this, but the basic idea
// is outlined below.
//
// Likely we would put these in a particular namespace so that the functions
// can have simple names, but then also not clash. Something like
// interface::tree::?
//
// This would make it simple to specialize for same nodes for S and T, or
// different nodes. The user just needs to specialize for the ones needed.
// The evaluator will still have two node types.

template <typename TreeNode>
int n_children(TreeNode *node);

template <typename TreeNode>
TreeNode *child(TreeNode *node, int i);

template <typename TreeNode>
int locality(TreeNode *node);

// Then the user provides specializations for using their tree.
template <>
SomeUserNode *child(int i) {
  return SomeUserNode->child[i];
}