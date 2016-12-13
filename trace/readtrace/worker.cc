#include "worker.h"

#include <algorithm>


namespace trace {


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
}


void Worker::add_file(File &stream) {
  auto endnow{events_.size()};

  do {
    events_.push_back(stream.next());
  } while (events_[events_.size() - 1] != nullptr);
  events_.pop_back();

  std::inplace_merge(events_.begin(), events_.begin() + endnow, events_.end(),
                     event_compare_uptr);
}


} // trace