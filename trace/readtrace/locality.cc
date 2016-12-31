#include "locality.h"

#include <exception>
#include <limits>


namespace traceutils {


Locality::Locality(File &stream)
    : locality_{stream.locality()}, workers_{}, locked_{false} {
  Worker to_add{stream};
  workers_[stream.worker()] = std::move(to_add);
}


std::map<int, int> Locality::worker_map() const {
  std::map<int, int> retval{};
  for (auto i = workers_.begin(); i != workers_.end(); ++i) {
    retval[i->first] = 0;
  }
  return retval;
}


size_t Locality::num_events() const {
  size_t retval{};
  for (auto i = workers_.begin(); i != workers_.end(); ++i) {
    retval += i->second.num_events();
  }
  return retval;
}


void Locality::add_file(File &stream) {
  if (locked_) {
    throw std::invalid_argument("Cannot add file to locked object");
  }
  if (stream.locality() != locality_) {
    throw std::invalid_argument("File does not match locality");
  }

  auto wkr = workers_.find(stream.worker());
  if (wkr != workers_.end()) {
    wkr->second.add_file(stream);
  } else {
    Worker to_add{stream};
    workers_[stream.worker()] = std::move(to_add);
  }
}


void Locality::zero_reset() {
  // First look over each worker - if there is one, it will be the first
  // event
  uint64_t zero{0};
  uint64_t negative_one = -1;
  for (auto i = workers_.begin(); i != workers_.end(); ++i) {
    zero = i->second.get_zeroref();
    if (zero != negative_one) {
      break;
    }
  }

fprintf(stdout, "DEBUG: locality %d has zero = %lu\n", locality_, zero);

  // Then if found, go ahead and offset all events in all workers
  // We also skip this work if zero is equal to 0.
  if (zero != negative_one && zero != 0) {
    for (auto i = workers_.begin(); i != workers_.end(); ++i) {
      i->second.reset_zero(zero);
    }
  }
}


void Locality::finalize(uint64_t min, uint64_t max) {
  if (locked_) return;

  locked_ = true;
  for (auto i = workers_.begin(); i != workers_.end(); ++i) {
    i->second.finalize(min, max);
  }
}


uint64_t Locality::max_ns() const {
  uint64_t retval{std::numeric_limits<uint64_t>::min()};
  for (auto iter = workers_.begin(); iter != workers_.end(); ++iter) {
    retval = std::max(retval, iter->second.max_ns());
  }
  return retval;
}


uint64_t Locality::min_ns() const {
  uint64_t retval{std::numeric_limits<uint64_t>::max()};
  for (auto iter = workers_.begin(); iter != workers_.end(); ++iter) {
    retval = std::min(retval, iter->second.min_ns());
  }
  return retval;
}


span_t Locality::span(uint64_t start, uint64_t end) const {
  if (!locked_) {
    throw std::invalid_argument("Cannot span unless locked");
  }

  span_t retval{};
  for (auto iter = workers_.begin(); iter != workers_.end(); ++iter) {
    retval[iter->first] = iter->second.range(start, end);
  }
  return retval;
}


} // traceutils
