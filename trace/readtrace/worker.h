#ifndef __TRACEUTILS_WORKER_H__
#define __TRACEUTILS_WORKER_H__


#include <memory>
#include <vector>

#include "traceevent.h"
#include "tracefile.h"


namespace traceutils {


/*

So this is going to be a container for events. This will basically want to
accept a file object and absorb the events in that file. The point will be
that the events get sorted into a hierarchical structure (a tree probably) by
their timestamp. This will make eventual use and searches easier.

One question is how to deal with the range of events. One observation is that
the trace files are more or less chronological. Or should we instead just use
a big damn vector, and apply over the top a tree structure to make it easier
to search? I think that is probably a good one.

*/


// This controls the refinement of the tree.
constexpr size_t kTreePartitionLimit = 40;


class Worker {
 public:
  Worker(File &stream);
  Worker(int id = -1, int loc = -1) noexcept 
      : id_{id}, loc_{loc}, events_{}, locked_{false}, root_{nullptr} { }
  Worker(Worker &&other) = default;

  // Make sure we get the default move
  Worker &operator=(Worker &&other) = default;

  // Simple queries
  int id() const {return id_;}
  int locality() const {return loc_;}
  size_t num_events() const {return events_.size();}

  // This could throw an exception std::runtime_error if the File object
  // throws one. this can throw std::invalid_argument if the stream's
  // worker ID does not match, or if the locality does not match.
  void add_file(File &stream);

  // This finalizes input and builds the tree
  void finalize(uint64_t min, uint64_t max);

  // Query the data in various ways
  // TODO: This is where we shall have to add the ability to search for
  // specific windows and things.
  uint64_t max_ns() const;
  uint64_t min_ns() const;

 private:
  struct node_t {
    using iter_t = std::vector<std::unique_ptr<Event>>::const_iterator;

    // These are the child nodes
    std::unique_ptr<node_t> left;
    std::unique_ptr<node_t> right;
    // The range in times
    uint64_t low;
    uint64_t high;
    // iterators into the events
    iter_t first;
    iter_t last;

    // Constructor
    node_t(uint64_t l, uint64_t h, iter_t f, iter_t l);
    // middle of the range in timestamps
    uint64_t middle() const;
    // partition a node
    void partition();
  };

  void build_tree(uint64_t min, uint64_t max);

  int id_;  // Which worker is this (numbered per locality)
  int loc_; // Which locality for which this was a worker
  std::vector<std::unique_ptr<Event>> events_;
  bool locked_;
  std::unique_ptr<node_t> root_;
};


} // traceutils


#endif // __TRACEUTILS_WORKER_H__