#ifndef __TRACEUTILS_SEGMENT_H__
#define __TRACEUTILS_SEGMENT_H__


#include <map>
#include <memory>
#include <vector>

#include "traceevent.h"


namespace traceutils {


class Segment {
 public:
  Segment()
      : id_{-1},
        begin_{-1},
        end_{-1},
        has_begin_{false},
        has_end_{false} { }

  bool claim(const Event *evt, int eid) {
    if (full()) {
      fprintf(stderr, "Segment is already full\n");
      return false;
    }

    if (!has_begin_) {
      id_ = evt->segment_type();
      if (!evt->start()) {
        fprintf(stderr, "ERROR: not a beginning of segment\n");
        return false;
      }
      begin_ = eid;
      has_begin_ = true;
    } else if (!has_end_) {
      if (evt->segment_type() != id_) {
        fprintf(stderr, "ERROR: segment type does not match\n");
        return false;
      }
      if (!evt->end()) {
        fprintf(stderr, "ERROR: not an end of segment\n");
        return false;
      }
      end_ = eid;
      has_end_ = true;
    }

    return true;
  }

  bool full() const {return has_begin_ && has_end_;}

  int segment() const {return id_;}
  int begin() const {return begin_;}
  int end() const {return end_;}

 private:
  int id_;
  int begin_;
  int end_;
  bool has_begin_;
  bool has_end_;
};


class SegmentCollator {
 public:
  SegmentCollator() : partial_{}, full_{} { }

  bool collate(const Event *evt, int eid) {
    int sid = evt->segment_type();
    if (sid == 0) return true;

    // add to map as needed
    if (partial_.find(sid) == partial_.end()) {
      partial_[sid] = std::unique_ptr<Segment>{nullptr};
    }

    // create a Segment as needed
    if (partial_[sid] == nullptr) {
      partial_[sid] = std::unique_ptr<Segment>{new Segment{}};
    }

    // claim the event
    if (!(partial_[sid])->claim(evt, eid)) {
      return false;
    }

    if (partial_[sid]->full()) {
      full_.push_back(std::move(partial_[sid]));
      partial_[sid] = nullptr;
    }

    return true;
  }

  // add segments from one to the other
  void combine(SegmentCollator &other) {
    for (size_t i = 0; i < other.full_.size(); ++i) {
      this->full_.push_back(std::move(other.full_[i]));
    }
  }

  size_t SegmentCount() const {return full_.size();}
  const Segment *GetSegment(size_t index) const {return full_[index].get();}

 private:
  std::map<int, std::unique_ptr<Segment>> partial_;
  std::vector<std::unique_ptr<Segment>> full_;
};


} // traceutils


#endif // __TRACEUTILS_SEGMENT_H__