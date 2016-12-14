#ifndef __TRACEUTILS_NETWORK_EVENTS_H__
#define __TRACEUTILS_NETWORK_EVENTS_H__


#include "traceevent.h"


namespace traceutils {

namespace network {


// The event types currently supported
enum class EventType {
  kProgressBegin,
  kProgressEnd,
  kUnknown
};


// Create a prototype for the file object
std::unique_ptr<Event> prototype_from_name(const std::string &name);


// The ProgressBegin event
class ProgressBegin : public Event {
 public:
  ProgressBegin(uint64_t stamp = 0) noexcept : Event{stamp} { }

  const std::string &event_class() const override;
  const std::string &event_type() const override;
  int num_fields() const noexcept override {return 0;}
  uint64_t field(int i) const noexcept override {return 0;}
  bool start() const noexcept override {return true;}
  bool end() const noexcept override {return false;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const override;
};


// The ProgressEnd event
class ProgressEnd : public Event {
 public:
  ProgressEnd(uint64_t stamp = 0) noexcept : Event{stamp} { }

  const std::string &event_class() const override;
  const std::string &event_type() const override;
  int num_fields() const noexcept override {return 0;}
  uint64_t field(int i) const noexcept override {return 0;}
  bool start() const noexcept override {return false;}
  bool end() const noexcept override {return true;}
  std::unique_ptr<Event> read_from_file(FILE *fd) const override;
};


} // traceutils::network

} // traceutils


#endif // __TRACEUTILS_NETWORK_EVENTS_H__