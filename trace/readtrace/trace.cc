#include "trace.h"


#include <limits>


namespace traceutils {


int Trace::num_workers(int loc) const {
  auto iter = locs_.find(loc);
  if (iter != locs_.end()) {
    return iter->second.num_workers();
  } else {
    return -1;
  }
}


int Trace::max_worker(int loc) const {
  auto iter = locs_.find(loc);
  if (iter != locs_.end()) {
    return iter->second.max_worker();
  } else {
    return -1;
  }
}

size_t Trace::num_events(int loc) const {
  auto iter = locs_.find(loc);
  if (iter != locs_.end()) {
    return iter->second.num_events();
  } else {
    return -1;
  }
}


size_t Trace::num_events() const {
  size_t retval{0};
  for (auto i = locs_.begin(); i != locs_.end(); ++i) {
    retval += i->second.num_events();
  }
  return retval;
}


void Trace::add_file(File &stream) {
  if (locked_) {
    throw std::invalid_argument("Cannot add file to locked object");
  }

  auto iter = locs_.find(stream.locality());
  if (iter != locs_.end()) {
    iter->second.add_file(stream);
  } else {
    Locality to_add{stream};
    locs_[stream.locality()] = std::move(to_add);
  }
}


void Trace::finalize() {
  if (locked_) return;

  locked_ = true;
  uint64_t begin = min_ns();
  uint64_t end = max_ns() + 1;
  for (auto i = locs_.begin(); i != locs_.end(); ++i) {
    i->second.finalize(begin, end);
  }
}


uint64_t Trace::max_ns() const {
  uint64_t retval{std::numeric_limits<uint64_t>::min()};
  for (auto iter = locs_.begin(); iter != locs_.end(); ++iter) {
    retval = std::max(retval, iter->second.max_ns());
  }
  return retval;
}


uint64_t Trace::min_ns() const {
  uint64_t retval{std::numeric_limits<uint64_t>::max()};
  for (auto iter = locs_.begin(); iter != locs_.end(); ++iter) {
    retval = std::min(retval, iter->second.min_ns());
  }
  return retval;
}


window_t Trace::window(uint64_t start, uint64_t end) const {
  if (!locked_) {
    throw std::invalid_argument("Cannot window unless locked");
  }

  window_t retval{};
  for (auto iter = locs_.begin(); iter != locs_.end(); ++iter) {
    retval[iter->first] = iter->second.span(start, end);
  }
  return retval;
}


} //traceutils
