#ifndef __TRACEUTILS_TRACE_H__
#define __TRACEUTILS_TRACE_H__


#include <map>

#include "locality.h"
#include "tracefile.h"


namespace traceutils {


using window_t = std::map<int, span_t>;


class Trace {
 public:
  Trace() : locked_{false}, locs_{} { }

  int localities() const {return locs_.size();}
  int max_locality() const {return locs_.rbegin()->first;}

  int num_workers(int loc) const;
  int max_worker(int loc) const;

  size_t num_events(int loc) const;
  size_t num_events() const;

  void add_file(File &stream);

  // This finalizes additions and creates the trees
  void finalize();

  // Get the min and max record time of any locality
  uint64_t max_ns() const;
  uint64_t min_ns() const;

  // Get a window of events
  window_t window(uint64_t start, uint64_t end) const;

 private:
  bool locked_;
  std::map<int, Locality> locs_;
};


} // traceutils


#endif // __TRACEUTILS_TRACE_H__