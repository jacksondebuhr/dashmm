#include "trace.h"


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
    throw // TODO: something appropriate
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
  uint64_t end = max_ns();
  for (auto i = locs_.begin(); i != locs_.end(); ++i) {
    i->finalize(begin, end);
  }
}


uint64_t Trace::max_ns() const {
  // TODO: check this limits thing
  uint64_t retval{std::numeric_limits<uint64_t>::min};
  for (auto iter = locs_.begin(); iter != locs_.end(); ++iter) {
    retval = std::max(retval, iter->second.max_ns());
  }
  return retval;
}


uint64_t Trace::min_ns() const {
  // TODO: check this limits thing
  uint64_t retval{std::numeric_limits<uint64_t>::max};
  for (auto iter = locs_.begin(); iter != locs_.end(); ++iter) {
    retval = std::min(retval, iter->second.min_ns());
  }
  return retval;
}


} //traceutils
