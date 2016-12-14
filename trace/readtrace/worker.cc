#include "worker.h"

#include <cassert>

#include <algorithm>
#include <exception>


namespace traceutils {


Worker::Worker(File &stream)
    : id_{stream.worker()}, loc_{stream.locality()}, events_{},
      locked_{false} {
  std::vector<std::unique_ptr<Event>> trial{};
  do {
    trial.push_back(stream.next());
  } while (trial[trial.size() - 1] != nullptr);
  trial.pop_back();
  // If we have made it here, there was no problem reading the file, so
  // move it into this object's events_
  events_ = std::move(trial);
  // TODO: of course, if an exception is thrown, the File object will have
  // been modified. Do we care?

  // Here we release that resource now that we are done.
  stream.close();
}


void Worker::add_file(File &stream) {
  if (locked_) {
    throw std::invalid_argument("File cannot be added to locked object");
  }

  if (stream.worker() != id_ || stream.locality() != loc_) {
    throw std::invalid_argument("File does not match worker");
  }

  auto endnow{events_.size()};

  do {
    events_.push_back(stream.next());
  } while (events_[events_.size() - 1] != nullptr);
  events_.pop_back();

  // Release the file object.
  stream.close();

  std::inplace_merge(events_.begin(), events_.begin() + endnow, events_.end(),
                     event_compare_uptr);
}


void Worker::finalize(uint64_t min, uint64_t max) {
  if (locked_) return;
  locked_ = true;
}


// Because we keep the events sorted, the min and max times must be from the
// start and end of the vector of events.
uint64_t Worker::max_ns() const {
  return events_[events_.size() - 1]->stamp();
}

uint64_t Worker::min_ns() const {
  return events_[0]->stamp();
}


range_t Worker::range(uint64_t start, uint64_t end) const {
  if (!locked_) {
    throw std::invalid_argument("Cannot range unless locked");
  }
  if (start > end) {
    throw std::invalid_argument("inverted intervals are not allowed");
  }

  auto first = events_.begin();
  auto last = events_.end();
  auto low = events_[0]->stamp();
  auto high = events_[events_.size() - 1]->stamp();

  // Deal with two easy cases first
  if (end < low) {
    return make_pair(first, first);
  }
  if (start >= high) {
    return make_pair(last, last);
  }

  // Map range onto possible range (doing so excludes no events)
  start = std::max(start, low);
  end = std::min(end, high);

  // Find the range
  auto comp_s = [&start](const std::unique_ptr<Event> &a) {
    return a->stamp() < start;
  };
  iter_t one = std::partition_point(first, last, comp_s);

  auto comp_e = [&end](const std::unique_ptr<Event> &a) {
    return a->stamp() < end;
  };
  iter_t two = std::partition_point(first, last, comp_e);

  return make_pair(one, two);
}


} // traceutils
