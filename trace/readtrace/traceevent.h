#ifndef __TRACEUTILS_EVENT_H__
#define __TRACEUTILS_EVENT_H__


#include <cstdint>

#include <memory>
#include <string>


namespace traceutils {


// The classes of events that HPX-5 provides
enum class EventClass {
  kParcel,
  kNetwork,
  kSched,
  kLCO,
  kProcess,
  kMemory,
  kTrace,
  kGAS,
  kCollective,
  kUnknown
};


EventClass event_class_from_name(const std::string &name);


// This is going to be an abstract base class from which the other events
// will derive. The one thing that is sure about all events is that it will
// have the uint64_t timestamp. So that will be part of the base class.
class Event {
 public:
  Event(uint64_t stamp = 0) noexcept : stamp_{stamp} { }
  virtual ~Event() { }


  uint64_t stamp() const noexcept {return stamp_;}
  void set_stamp(uint64_t s) noexcept {stamp_ = s;}

  // What class of event is this
  virtual const std::string &event_class() const = 0;

  // What type of event is this
  virtual const std::string &event_type() const = 0;

  // How many extra fields are present
  virtual int num_fields() const noexcept = 0;

  // Retrieve the extra fields
  virtual uint64_t field(int i) const noexcept = 0;

  // A factory method that will produce a new event from a given file
  //
  // This will emit an exception at the end of file during a incomplet read
  // This will emit an exception when there is some IO problem
  //
  // If the end of file is reached at the start of a read, a null pointer is
  // returned.
  virtual std::unique_ptr<Event> read_from_file(FILE *fd) const = 0;

 private:
  // All events have a timestam
  uint64_t stamp_;
};


// Comparison for unique pointer case
bool event_compare_uptr(const std::unique_ptr<Event> &a,
                        const std::unique_ptr<Event> &b);


} //traceutils


#endif // __TRACEUTILS_EVENT_H__