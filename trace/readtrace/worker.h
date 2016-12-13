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
  // TODO: A constructor taking a File object?
  Worker(File &stream);
  Worker(int id, int loc) noexcept : id_{id}, loc_{loc}, events_{} { }

  // Simple queries
  int id() const {return id_;}
  int locality() const {return loc_;}
  size_t num_events() const {return events_.size();}

  // This could throw an exception std::runtime_error if the File object
  // throws one.
  void add_file(File &stream);

  // Query the data in various ways

private:
  int id_;  // Which worker is this (numbered per locality)
  int loc_; // Which locality for which this was a worker
  std::vector<std::unique_ptr<Event>> events_;
};


} // trace


#endif // __TRACE_WORKER_H__