// Some C++ like pseudo code for the tree separation

// This set of functions deals with discovering the DAG component confined to
// the source tree
//
/*

So the DualTree which in the current scheme is mostly meta data would be
absorbed into the Evaluator object, which is the most sensible place for these
functions to go.

Then, DualTree would really have no part. instead, anything from a tree, which
would ultimately be in the purview of a tree provider, would need to be
made available through one of the template interface functions.

In the following, things that are written as node->*** will be made into an
element of that template interface.

The DAG will be part of the meta object, and will have both collective access
to the nodes, and also will have maps from tree node addresses into the
DAG nodes.

*/
SourceApplyMethod(Evaluator *meta, SourceNode *node, LCO done) {
  int n_children = interface::tree::n_children(node)
  if (n_children == 0) {
    // Leaves generate
    DAGNode *P, *N, *I;
    (P, N, I) = Evaluator::method_t::generate(node, meta->domain());
    int local = interface::tree::locality(node);
    P->locality = local;
    N->locality = local;
    Evaluator::method_t::distropolicy_t::assign_for_source(P, N, I, local);
    meta->dag.add(P, N, I, node);
    done.set();
  } else {
    // Internal nodes aggregate, after the children are done
    LCO cdone = new LCO::And(n_children);
    for (int i = 0; i < 8; ++i) {
      SourceNode *child = interface::tree::child(node, i);
      if (child) {
        call SourceApplyMethod(meta, child, cdone);
      }
    }

    when (cdone) {
      call SourceApplyMethodChildDone(meta, node, done);
    } then {
      delete cdone;
    }
  }

  return;
}

SourceApplyMethodChildDone(DualTree *meta, SourceNode *node, LCO done) {
  DAGNode *N, *I;
  (N, I) = Evaluator::method_t::aggregate(node, meta->domain());
  Evaluator::method_t::distropolicy_t::assign_for_source(
        nullptr, N, I, interface::tree::locality(node));
  meta->dag.add(nullptr, N, I, node);
  done.set();
  return;
}

