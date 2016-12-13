#include "worker.h"

#include <algorithm>
#include <exception>


namespace traceutils {


Worker::Worker(File &stream)
    : id_{stream.worker()}, loc_{stream.locality()}, events_{} {
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


} // traceutils