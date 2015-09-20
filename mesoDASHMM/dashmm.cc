#include "dashmm.h"

#include <algorithm>
#include <cassert>
#include <memory>
#include <map>
#include <vector>

#include "BH_builtins.h"
#include "BH_expansion.h"
#include "BH_method.h"
#include "domaingeometry.h"
#include "expansion.h"
#include "FMM_builtins.h"
#include "FMM_expansion.h"
#include "FMM_method.h"
#include "method.h"
#include "node.h"


namespace dashmm {

//NOTE: The only thing that should appear in this file should be the
// interfaces. Extra utility functions should appear elsewhere to keep this
// file clean, and easy to maintain. (This policy subject to change).


void init() {
  bh_builtins_init();
  fmm_builtins_init();
}


void evaluate(std::vector<Source> &sources,
              std::vector<Target> &targets,
              double limit,
              Method *method,
              Expansion *expand) {
  if (method == nullptr) {
    return;
  }

  //Probably factor this, eh?
  //NOTE: we are using the same root for both trees. This is probably not too
  // much waste overall. So I will not worry too much about doing
  //something more sophisticated here.

  double min[3]{1.0e34, 1.0e34, 1.0e34};
  double max[3]{-1.0e34, -1.0e34, -1.0e34};
  for (auto i = sources.begin(); i != sources.end(); ++i) {
    if (i->x() < min[0]) min[0] = i->x();
    if (i->y() < min[1]) min[1] = i->y();
    if (i->z() < min[2]) min[2] = i->z();
    if (i->x() > max[0]) max[0] = i->x();
    if (i->y() > max[1]) max[1] = i->y();
    if (i->z() > max[2]) max[2] = i->z();
  }
  for (auto i = targets.begin(); i != targets.end(); ++i) {
    if (i->x() < min[0]) min[0] = i->x();
    if (i->y() < min[1]) min[1] = i->y();
    if (i->z() < min[2]) min[2] = i->z();
    if (i->x() > max[0]) max[0] = i->x();
    if (i->y() > max[1]) max[1] = i->y();
    if (i->z() > max[2]) max[2] = i->z();
    i->set_phi(0.0);
  }
  DomainGeometry geo{Point{min[0], min[1], min[2]},
      Point{max[0], max[1], max[2]}};

  SourceNode *source_root = new SourceNode(geo, 0, 0, 0, 0, method, nullptr);
  TargetNode *target_root = new TargetNode(geo, 0, 0, 0, 0, method, nullptr);

  source_root->partition(sources.begin(), sources.end(), limit, expand);

  std::vector<SourceNode *> consider(1);
  consider[0] = source_root;
  target_root->partition(targets.begin(), targets.end(), limit, expand, 0,
                         consider);


  //target_root->descend(expand, 0, consider);

  delete source_root;
  delete target_root;
}


void finalize() {
  bh_builtins_fini();
  fmm_builtins_fini();
}


} //namespace dashmm
