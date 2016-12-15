#include "trace.h"

#include <cassert>

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


int Trace::total_workers() const {
  int retval{0};
  for (auto iter = locs_.begin(); iter != locs_.end(); ++iter) {
    retval += iter->second.num_workers();
  }
  return retval;
}


std::map<int, std::map<int, int>> Trace::worker_map() const {
  std::map<int, std::map<int, int>> retval{};
  for (auto i = locs_.begin(); i != locs_.end(); ++i) {
    retval[i->first] = i->second.worker_map();
  }
  return retval;
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


namespace {
  double range_coverage(const range_t &range, int type,
                        uint64_t start, uint64_t end) {
    // Shortcut for no entries
    if (range.second - range.first == 0) return 0.0;

    uint64_t total{end - start};
    uint64_t cover{0};
    int64_t previous_start = start;
    int64_t previous_end = start - 1;
    for (auto i = range.first; i != range.second; ++i) {
      if ((*i)->segment_type() == type) {
        if ((*i)->start()) {
          // make sure the segments are alternating start/stop
          previous_start = (*i)->stamp();
          assert(previous_end < previous_start);
        } else {
          // make sure if we are not a start we are an end
          assert((*i)->end());
          // make sure the segments are alternating start and end
          previous_end = (*i)->stamp();
          assert(previous_end >= previous_start);
          cover += previous_end - previous_start;
        }
      }
    }
    // Deal with a start but no end here
    if (previous_end < previous_start) {
      cover += end - previous_start;
    }
    return ((double)cover) / total;
  }

  std::map<int, double> span_coverage(const span_t &span, int type,
                                      uint64_t start, uint64_t end) {
    std::map<int, double> retval{};
    for (auto i = span.begin(); i != span.end(); ++i) {
      retval[i->first] = range_coverage(i->second, type, start, end);
    }
    return retval;
  }
} // anonymous


cover_t coverage_of_segment_type(const window_t &window, int type,
                                 uint64_t start, uint64_t end) {
  cover_t retval{};
  for (auto i = window.begin(); i != window.end(); ++i) {
    retval[i->first] = span_coverage(i->second, type, start, end);
  }
  return retval;
}


} //traceutils
