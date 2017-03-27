#ifndef __TRACEUTILS_UTILIZATION_H__
#define __TRACEUTILS_UTILIZATION_H__


#include <string>
#include <vector>

#include "trace.h"


namespace traceutils {


class Utilization {
 public:
  Utilization(const Trace &trace) : trace_{trace} { }

  // Create a utilization plot data
  //
  // fname - outfile name
  // t0 - initial timestamp
  // t1 - final timestamp
  // samples - the number of bins in time
  // segments - the segments to sum over; this will create one column for each
  //            segment, as well as one for the sum of the input segments
  void operator()(const std::string &fname, uint64_t t0, uint64_t t1,
                  int samples, const std::vector<int> segments);

 private:
  const Trace &trace_;
};


} // traceutils


#endif // __TRACEUTILS_UTILIZATION_H__