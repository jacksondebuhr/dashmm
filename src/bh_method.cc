#include "include/bh_method.h"

//C/C++

// other DASHMM


namespace dashmm {


//


Method *bh_method(double theta) {
  Method *retval = new BHMethod{theta};
  return retval;
}


} // namespace dashmm
