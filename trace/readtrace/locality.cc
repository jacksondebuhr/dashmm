#include "locality.h"

#include <exception>
#include <limits>


namespace traceutils {


Locality::Locality(File &stream) 
    : locality_{stream.locality()}, workers_{}, locked_{false} {
  Worker to_add{stream};
  workers_[stream.worker()] = std::move(to_add);
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
    throw // TODO: something appropriate
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


void Locality::finalize(uint64_t min, uint64_t max) {
  if (locked_) return;

  locked_ = true;
  for (auto i = workers_.begin(); i != workers_.end(); ++i) {
    i->second.finalize(min, max);
  }
}


uint64_t Locality::max_ns() const {
  uint64_t retval{std::numeric_limits<uint64_t>::min};
  for (auto iter = workers_.begin(); iter != workers_.end(); ++iter) {
    retval = std::max(retval, iter->second.max_ns());
  }
  return retval;
}


uint64_t Locality::min_ns() const {
  uint64_t retval{std::numeric_limits<uint64_t>::max};
  for (auto iter = workers_.begin(); iter != workers_.end(); ++iter) {
    retval = std::min(retval, iter->second.min_ns());
  }
  return retval;
}


} // traceutils
