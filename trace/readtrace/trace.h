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
  int total_workers() const;
  // TODO
  // Same idea. Make it functional, taking a template parameter which sets
  // the value.
  std::map<int, std::map<int, int>> worker_map() const;

  size_t num_events(int loc) const;
  size_t num_events() const;

  void add_file(File &stream);

  // For each locality, find the zeroref event and if it exists, shift each
  // event to start from the given zero.
  void find_and_reset_zero();

  // This finalizes additions
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


// The first is locality, the second is worker
using cover_t = std::map<int, std::map<int, double>>;

cover_t coverage_of_segment_type(const window_t &window, int type,
                                 uint64_t start, uint64_t end);


} // traceutils


#endif // __TRACEUTILS_TRACE_H__