#include "networkevents.h"

#include <exception>


namespace trace {

namespace network {

namespace {
  // The class of events implemented here
  static std::string kEventClass{"network"};

  // The various event types
  static std::string kProgressBeginType{"progress begin"};
  static std::string kProgressEndType{"progress end"};

  // Convert name to type - this may not be needed elsewhere, and so we
  // may drop this eventually.
  EventType type_from_name(const std::string &name) {
    if (name == std::string("PROGRESS_BEGIN")) {
      return EventType::kProgressBegin;
    } else if (name == std::string("PROGRESS_END")) {
      return EventType::kProgressEnd;
    } else {
      return EventType::kUnknown;
    }
  }

} // anonymous


std::unique_ptr<Event> prototype_from_name(const std::string &name) {
  auto type = type_from_name(name);

  switch (type) {
    case EventType::kProgressBegin:
      return std::unique_ptr<Event>{new ProgressBegin{}};
      break;
    case EventType::kProgressEnd:
      return std::unique_ptr<Event>{new ProgressEnd{}};
      break;
    case EventType::kUnknown:
      return std::unique_ptr<Event>{nullptr};
      break;
  }

  return std::unique_ptr<Event>{nullptr};
}


const std::string &ProgressBegin::event_class() const {
  return kEventClass;
}

const std::string &ProgressBegin::event_type() const {
  return kProgressBeginType;
}

std::unique_ptr<Event> ProgressBegin::read_from_file(FILE *fd) const {
  uint64_t rval{0};
  auto before = ftell(fd);
  if (1 != fread(&rval, sizeof(rval), 1, fd)) {
    auto after = ftell(fd);
    if (feof(fd)) {
      if (before == after) {
        return std::unique_ptr<Event>{nullptr};
      } else {
        throw std::runtime_error("incomplete read, suggesting format error");
      }
    } else if (ferror(fd)) {
      // TODO make this a bit more meaningful
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new ProgressBegin{rval}};
}




const std::string &ProgressEnd::event_class() const {
  return kEventClass;
}

const std::string &ProgressEnd::event_type() const {
  return kProgressEndType;
}

std::unique_ptr<Event> ProgressEnd::read_from_file(FILE *fd) const {
  uint64_t rval{0};
  auto before = ftell(fd);
  if (1 != fread(&rval, sizeof(rval), 1, fd)) {
    auto after = ftell(fd);
    if (feof(fd)) {
      if (before == after) {
        return std::unique_ptr<Event>{nullptr};
      } else {
        throw std::runtime_error("incomplete read, suggesting format error");
      }
    } else if (ferror(fd)) {
      // TODO make this a bit more meaningful
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new ProgressEnd{rval}};
}



} // trace::network

} // trace