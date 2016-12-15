#ifndef __TRACEUTILS_LOCALITY_H__
#define __TRACEUTILS_LOCALITY_H__


#include <map>

#include "tracefile.h"
#include "worker.h"


namespace traceutils {


using span_t = std::map<int, range_t>;


class Locality {
 public:
  Locality(File &stream);
  Locality(int loc = -1) : locality_{loc}, workers_{}, locked_{false} { }
  Locality(Locality &&other) = default;

  // Make sure default move operator exists
  Locality &operator=(Locality &&other) = default;


  int locality() const noexcept {return locality_;}
  // The number of workers currently.
  int num_workers() const {return workers_.size();}
  // The maximum id of any workers.
  int max_worker() const {return workers_.rbegin()->first;}
  // Get a map from worker id to integers. This is initialized with zeros
  // Maybe make this functional by taking a template parameter for a callable
  // type. That could be the thing that sets the value.
  std::map<int, int> worker_map() const;

  size_t num_events() const;

  // This could throw std::runtime_error or std::invalid_argument
  void add_file(File &stream);

  // This finalizes input and builds the trees
  void finalize(uint64_t min, uint64_t max);

  // Min and max times in this Locality
  uint64_t max_ns() const;
  uint64_t min_ns() const;

  // Grab a span
  span_t span(uint64_t start, uint64_t end) const;

 private:
  int locality_;
  std::map<int, Worker> workers_;
  bool locked_;
};


} // traceutils


#endif // __TRACEUTILS_LOCALITY_H__