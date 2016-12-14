#ifndef __TRACEUTILS_TRACE_H__
#define __TRACEUTILS_TRACE_H__


#include <map>

#include "locality.h"
#include "tracefile.h"


namespace traceutils {


class Trace {
 public:
  Trace() : locked_{false}, locs_{} { }

  int localities() const {return locs_.size();}
  int max_locality() const {return locs_.rbegin()->first;}

  int num_workers(int loc) const;
  int max_worker(int loc) const;

  // TODO do we want both here?
  size_t num_events(int loc) const;
  size_t num_events() const;

  void add_file(File &stream);

  // This finalizes additions and creates the trees
  void finalize();

  // Get the min and max record time of any locality
  uint64_t max_ns() const;
  uint64_t min_ns() const;
 private:
  bool locked_;
  std::map<int, Locality> locs_;
};


} // traceutils


#endif // __TRACEUTILS_TRACE_H__