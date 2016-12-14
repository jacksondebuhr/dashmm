#include "worker.h"

#include <algorithm>
#include <exception>


namespace traceutils {


Worker::Worker(File &stream)
    : id_{stream.worker()}, loc_{stream.locality()}, events_{}, 
      locked_{false}, root_{nullptr} {
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
    throw // TODO: pick an appropriate exception
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
  build_tree(min, max);
}


// Because we keep the events sorted, the min and max times must be from the
// start and end of the vector of events.
uint64_t Worker::max_ns() const {
  return events_[events.size() - 1]->stamp();
}

uint64_t Worker::min_ns() const {
  return events_[0]->stamp();
}


void Worker::build_tree(uint64_t min, uint64_t max) {
  root_ = std::unique_ptr<node_t>{
                  new node_t{min, max, events_.begin(), events_.end()}};
  root_->partition();
}


Worker::node_t::node_t(uint64_t l, uint64_t h, iter_t f, iter_t l) 
    : left{nullptr}, right{nullptr}, low{l}, high{h}, first{f}, last{l} { }


uint64_t Worker::node_t::middle() const {
  return (low + high) / 2;
}


void Worker::node_t::partition() {
  // Are we done?
  if (last - first <= kTreePartitionLimit) {
    return;
  }

  // if not, find split points
  uint64_t midpoint = middle();
  auto comp = [&middle](const std::unique_ptr<Event> &a) {
    return a->stamp() < midpoint;
  }
  iter_t split = std::partition_point(first, last, comp);

  // create new nodes
  if (split - first) {
    left = std::unique_ptr<node_t>{new node_t{low, middle, first, split}};
    left->partition();
  }
  if (last - split) {
    right = std::unique_ptr<node_t>{new node_t{middle, high, split, last}};
    right->partition();
  }
}


} // traceutils