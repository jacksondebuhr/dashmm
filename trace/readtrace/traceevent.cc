#include "traceevent.h"


namespace trace {


EventClass event_class_from_name(const std::string &name) {
  if (name == "PARCEL") {
    return EventClass::kParcel;
  } else if (name == "NETWORK") {
    return EventClass::kNetwork;
  } else if (name == "SCHED") {
    return EventClass::kSched;
  } else if (name == "LCO") {
    return EventClass::kLCO;
  } else if (name == "PROCESS") {
    return EventClass::kProcess;
  } else if (name == "MEMORY") {
    return EventClass::kMemory;
  } else if (name == "TRACE") {
    return EventClass::kTrace;
  } else if (name == "GAS") {
    return EventClass::kGAS;
  } else if (name == "COLLECTIVE") {
    return EventClass::kCollective;
  } else {
    return EventClass::kUnknown;
  }
}


bool event_compare_uptr(const std::unique_ptr<Event> &a,
                        const std::unique_ptr<Event> &b) {
  return a->stamp() < b->stamp();
}


} // trace
