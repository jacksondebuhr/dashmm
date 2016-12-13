#ifndef __TRACE_WORKER_H__
#define __TRACE_WORKER_H__


#include <memory>
#include <vector>

#include "traceevent.h"
#include "tracefile.h"


namespace trace {


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


class Worker {
 public:
  Worker(File &stream);
  Worker(int id = -1, int loc = -1) noexcept : id_{id}, loc_{loc}, events_{} { }

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

  // Query the data in various ways
  // TODO: This is where we shall have to add the ability to search for
  // specific windows and things.
 private:
  int id_;  // Which worker is this (numbered per locality)
  int loc_; // Which locality for which this was a worker
  std::vector<std::unique_ptr<Event>> events_;
};


} // trace


#endif // __TRACE_WORKER_H__