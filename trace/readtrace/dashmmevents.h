#ifndef __TRACEUTILS_DASHMM_EVENTS_H__
#define __TRACEUTILS_DASHMM_EVENTS_H__


#include "traceevent.h"


namespace traceutils {


// The segment identifiers for the events defined here
constexpr int kDASHMMStoT = 2;
constexpr int kDASHMMMtoT = 3;
constexpr int kDASHMMLtoT = 4;
constexpr int kDASHMMStoM = 5;
constexpr int kDASHMMStoL = 6;
constexpr int kDASHMMELCO = 7;
constexpr int kDASHMMMtoM = 8;
constexpr int kDASHMMMtoL = 9;
constexpr int kDASHMMLtoL = 10;
constexpr int kDASHMMMtoI = 11;
constexpr int kDASHMMItoI = 12;
constexpr int kDASHMMItoL = 13;


namespace dashmm {


enum class EventType {
  kStoTBegin,
  kStoTEnd,
  kMtoTBegin,
  kMtoTEnd,
  kLtoTBegin,
  kLtoTEnd,
  kStoMBegin,
  kStoMEnd,
  kStoLBegin,
  kStoLEnd,
  kELCOBegin,
  kELCOEnd,
  kMtoMBegin,
  kMtoMEnd,
  kMtoLBegin,
  kMtoLEnd,
  kLtoLBegin,
  kLtoLEnd,
  kMtoIBegin,
  kMtoIEnd,
  kItoIBegin,
  kItoIEnd,
  kItoLBegin,
  kItoLEnd,
  kZeroRef,
  kUnknown
};


std::unique_ptr<Event> prototype_from_name(const std::string &name);


class DASHMMEvent : public Event {
 public:
  DASHMMEvent(uint64_t stamp = 0) noexcept : Event{stamp} { }

  const std::string &event_class() const final;
  int num_fields() const noexcept final {return 0;}
  uint64_t field(int i) const noexcept final {return 0;}
};


// StoT
class StoTBegin : public DASHMMEvent {
 public:
  StoTBegin(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return true;}
  bool end() const noexcept final {return false;}
  int segment_type() const noexcept final {return kDASHMMStoT;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};

class StoTEnd : public DASHMMEvent {
 public:
  StoTEnd(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return false;}
  bool end() const noexcept final {return true;}
  int segment_type() const noexcept final {return kDASHMMStoT;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};


// MtoT
class MtoTBegin : public DASHMMEvent {
 public:
  MtoTBegin(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return true;}
  bool end() const noexcept final {return false;}
  int segment_type() const noexcept final {return kDASHMMMtoT;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};

class MtoTEnd : public DASHMMEvent {
 public:
  MtoTEnd(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return false;}
  bool end() const noexcept final {return true;}
  int segment_type() const noexcept final {return kDASHMMMtoT;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};


// LtoT
class LtoTBegin : public DASHMMEvent {
 public:
  LtoTBegin(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return true;}
  bool end() const noexcept final {return false;}
  int segment_type() const noexcept final {return kDASHMMLtoT;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};

class LtoTEnd : public DASHMMEvent {
 public:
  LtoTEnd(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return false;}
  bool end() const noexcept final {return true;}
  int segment_type() const noexcept final {return kDASHMMLtoT;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};


// StoM
class StoMBegin : public DASHMMEvent {
 public:
  StoMBegin(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return true;}
  bool end() const noexcept final {return false;}
  int segment_type() const noexcept final {return kDASHMMStoM;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};

class StoMEnd : public DASHMMEvent {
 public:
  StoMEnd(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return false;}
  bool end() const noexcept final {return true;}
  int segment_type() const noexcept final {return kDASHMMStoM;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};


// StoL
class StoLBegin : public DASHMMEvent {
 public:
  StoLBegin(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return true;}
  bool end() const noexcept final {return false;}
  int segment_type() const noexcept final {return kDASHMMStoL;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};

class StoLEnd : public DASHMMEvent {
 public:
  StoLEnd(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return false;}
  bool end() const noexcept final {return true;}
  int segment_type() const noexcept final {return kDASHMMStoL;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};


// ELCO
class ELCOBegin : public DASHMMEvent {
 public:
  ELCOBegin(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return true;}
  bool end() const noexcept final {return false;}
  int segment_type() const noexcept final {return kDASHMMELCO;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};

class ELCOEnd : public DASHMMEvent {
 public:
  ELCOEnd(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return false;}
  bool end() const noexcept final {return true;}
  int segment_type() const noexcept final {return kDASHMMELCO;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};


// MtoM
class MtoMBegin : public DASHMMEvent {
 public:
  MtoMBegin(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return true;}
  bool end() const noexcept final {return false;}
  int segment_type() const noexcept final {return kDASHMMMtoM;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};

class MtoMEnd : public DASHMMEvent {
 public:
  MtoMEnd(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return false;}
  bool end() const noexcept final {return true;}
  int segment_type() const noexcept final {return kDASHMMMtoM;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};


// MtoL
class MtoLBegin : public DASHMMEvent {
 public:
  MtoLBegin(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return true;}
  bool end() const noexcept final {return false;}
  int segment_type() const noexcept final {return kDASHMMMtoL;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};

class MtoLEnd : public DASHMMEvent {
 public:
  MtoLEnd(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return false;}
  bool end() const noexcept final {return true;}
  int segment_type() const noexcept final {return kDASHMMMtoL;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};


// LtoL
class LtoLBegin : public DASHMMEvent {
 public:
  LtoLBegin(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return true;}
  bool end() const noexcept final {return false;}
  int segment_type() const noexcept final {return kDASHMMLtoL;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};

class LtoLEnd : public DASHMMEvent {
 public:
  LtoLEnd(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return false;}
  bool end() const noexcept final {return true;}
  int segment_type() const noexcept final {return kDASHMMLtoL;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};


// MtoI
class MtoIBegin : public DASHMMEvent {
 public:
  MtoIBegin(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return true;}
  bool end() const noexcept final {return false;}
  int segment_type() const noexcept final {return kDASHMMMtoI;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};

class MtoIEnd : public DASHMMEvent {
 public:
  MtoIEnd(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return false;}
  bool end() const noexcept final {return true;}
  int segment_type() const noexcept final {return kDASHMMMtoI;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};


// ItoI
class ItoIBegin : public DASHMMEvent {
 public:
  ItoIBegin(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return true;}
  bool end() const noexcept final {return false;}
  int segment_type() const noexcept final {return kDASHMMItoI;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};

class ItoIEnd : public DASHMMEvent {
 public:
  ItoIEnd(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return false;}
  bool end() const noexcept final {return true;}
  int segment_type() const noexcept final {return kDASHMMItoI;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};


// ItoL
class ItoLBegin : public DASHMMEvent {
 public:
  ItoLBegin(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return true;}
  bool end() const noexcept final {return false;}
  int segment_type() const noexcept final {return kDASHMMItoL;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};

class ItoLEnd : public DASHMMEvent {
 public:
  ItoLEnd(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return false;}
  bool end() const noexcept final {return true;}
  int segment_type() const noexcept final {return kDASHMMItoL;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};


// ZeroRef
class ZeroRef : public DASHMMEvent {
 public:
  ZeroRef(uint64_t stamp = 0) noexcept : DASHMMEvent{stamp} { }

  const std::string &event_type() const final;
  bool start() const noexcept final {return false;}
  bool end() const noexcept final {return false;}
  int segment_type() const noexcept final {return 0;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const final;
};


} // traceutils::dashmm

} //traceutils


#endif // __TRACEUTILS_DASHMM_EVENTS_H__