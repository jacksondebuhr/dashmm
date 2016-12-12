#ifndef __TRACE_NETWORK_EVENTS_H__
#define __TRACE_NETWORK_EVENTS_H__


#include "traceevent.h"


namespace trace {

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
  std::unique_ptr<Event> read_from_file(FILE *fd) const override;
};


} // trace::network

} // trace


#endif // __TRACE_NETWORK_EVENTS_H__