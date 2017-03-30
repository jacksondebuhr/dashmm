#ifndef __TRACEUTILS_WORKER_H__
#define __TRACEUTILS_WORKER_H__


#include <memory>
#include <vector>

#include "traceevent.h"
#include "tracefile.h"


namespace traceutils {


// This controls the refinement of the tree.
constexpr int64_t kTreePartitionLimit = 10;

// Define these types
using iter_t = std::vector<std::unique_ptr<Event>>::const_iterator;
using range_t = std::pair<iter_t, iter_t>;


class Worker {
 public:
  Worker(File &stream);
  Worker(int id = -1, int loc = -1) noexcept
      : id_{id}, loc_{loc}, gid_{counter_++}, events_{}, locked_{false} { }
  Worker(Worker &&other) = default;

  // Make sure we get the default move
  Worker &operator=(Worker &&other) = default;

  // Simple queries
  int id() const {return id_;}
  int locality() const {return loc_;}
  int gid() const {return gid_;}
  size_t num_events() const {return events_.size();}

  // This could throw an exception std::runtime_error if the File object
  // throws one. this can throw std::invalid_argument if the stream's
  // worker ID does not match, or if the locality does not match.
  void add_file(File &stream);

  // This finalizes input and builds the tree
  void finalize(uint64_t min, uint64_t max);

  // Query the data in various ways
  uint64_t max_ns() const;
  uint64_t min_ns() const;

  // Get a range of events
  range_t range(uint64_t start, uint64_t end) const;

  // Get the zeroref if present; returns -1 is not present
  uint64_t get_zeroref() const;

  // Do something with all of the events
  template <typename Callable>
  void apply(Callable action) const {
    for (size_t i = 0; i < events_.size(); ++i) {
      action(events_[i].get());
    }
  }

  void reset_zero(uint64_t zero);

 private:
  int id_;  // Which worker is this (numbered per locality)
  int loc_; // Locality for which this was a worker
  int gid_;
  std::vector<std::unique_ptr<Event>> events_;
  bool locked_;
  static int counter_;
};


} // traceutils


#endif // __TRACEUTILS_WORKER_H__