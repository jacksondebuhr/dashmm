#include "dashmmevents.h"

#include <exception>


namespace traceutils {

namespace dashmm {

namespace {
  // The class of events implemented here
  static std::string kEventClass{"trace"};

  // The various event types
  static std::string kStoTBeginType{"dashmm stot begin"};
  static std::string kStoTEndType{"dashmm stot end"};
  static std::string kMtoTBeginType{"dashmm mtot begin"};
  static std::string kMtoTEndType{"dashmm mtot end"};
  static std::string kLtoTBeginType{"dashmm ltot begin"};
  static std::string kLtoTEndType{"dashmm ltot end"};
  static std::string kStoMBeginType{"dashmm stom begin"};
  static std::string kStoMEndType{"dashmm stom end"};
  static std::string kStoLBeginType{"dashmm stol begin"};
  static std::string kStoLEndType{"dashmm stl end"};
  static std::string kELCOBeginType{"dashmm elco begin"};
  static std::string kELCOEndType{"dashmm elco end"};
  static std::string kMtoMBeginType{"dashmm mtom begin"};
  static std::string kMtoMEndType{"dashmm mtom end"};
  static std::string kMtoLBeginType{"dashmm mtol begin"};
  static std::string kMtoLEndType{"dashmm mtol end"};
  static std::string kLtoLBeginType{"dashmm ltol begin"};
  static std::string kLtoLEndType{"dashmm ltol end"};
  static std::string kMtoIBeginType{"dashmm mtoi begin"};
  static std::string kMtoIEndType{"dashmm mtoi end"};
  static std::string kItoIBeginType{"dashmm itoi begin"};
  static std::string kItoIEndType{"dashmm itoi end"};
  static std::string kItoLBeginType{"dashmm itol begin"};
  static std::string kItoLEndType{"dashmm itol end"};
  static std::string kZeroRefType{"dashmm zero ref"};

  EventType type_from_name(const std::string &name) {
    if (name == std::string("DASHMM_STOT_BEGIN")) {
      return EventType::kStoTBegin;
    } else if (name == std::string("DASHMM_STOT_END")) {
      return EventType::kStoTEnd;
    } else if (name == std::string("DASHMM_MTOT_BEGIN")) {
      return EventType::kMtoTBegin;
    } else if (name == std::string("DASHMM_MTOT_END")) {
      return EventType::kMtoTEnd;
    } else if (name == std::string("DASHMM_LTOT_BEGIN")) {
      return EventType::kLtoTBegin;
    } else if (name == std::string("DASHMM_LTOT_END")) {
      return EventType::kLtoTEnd;
    } else if (name == std::string("DASHMM_STOM_BEGIN")) {
      return EventType::kStoMBegin;
    } else if (name == std::string("DASHMM_STOM_END")) {
      return EventType::kStoMEnd;
    } else if (name == std::string("DASHMM_STOL_BEGIN")) {
      return EventType::kStoLBegin;
    } else if (name == std::string("DASHMM_STOL_END")) {
      return EventType::kStoLEnd;
    } else if (name == std::string("DASHMM_ELCO_BEGIN")) {
      return EventType::kELCOBegin;
    } else if (name == std::string("DASHMM_ELCO_END")) {
      return EventType::kELCOEnd;
    } else if (name == std::string("DASHMM_MTOM_BEGIN")) {
      return EventType::kMtoMBegin;
    } else if (name == std::string("DASHMM_MTOM_END")) {
      return EventType::kMtoMEnd;
    } else if (name == std::string("DASHMM_MTOL_BEGIN")) {
      return EventType::kMtoLBegin;
    } else if (name == std::string("DASHMM_MTOL_END")) {
      return EventType::kMtoLEnd;
    } else if (name == std::string("DASHMM_LTOL_BEGIN")) {
      return EventType::kLtoLBegin;
    } else if (name == std::string("DASHMM_LTOL_END")) {
      return EventType::kLtoLEnd;
    } else if (name == std::string("DASHMM_MTOI_BEGIN")) {
      return EventType::kMtoIBegin;
    } else if (name == std::string("DASHMM_MTOI_END")) {
      return EventType::kMtoIEnd;
    } else if (name == std::string("DASHMM_ITOI_BEGIN")) {
      return EventType::kItoIBegin;
    } else if (name == std::string("DASHMM_ITOI_END")) {
      return EventType::kItoIEnd;
    } else if (name == std::string("DASHMM_ITOL_BEGIN")) {
      return EventType::kItoLBegin;
    } else if (name == std::string("DASHMM_ITOL_END")) {
      return EventType::kItoLEnd;
    } else if (name == std::string("DASHMM_ZEROREF")) {
      return EventType::kZeroRef;
    } else {
      return EventType::kUnknown;
    }
  }
} // anonymous


std::unique_ptr<Event> prototype_from_name(const std::string &name) {
  auto type = type_from_name(name);

  switch (type) {
    case EventType::kStoTBegin:
      return std::unique_ptr<Event>{new StoTBegin{}};
      break;
    case EventType::kStoTEnd:
      return std::unique_ptr<Event>{new StoTEnd{}};
      break;
    case EventType::kMtoTBegin:
      return std::unique_ptr<Event>{new MtoTBegin{}};
      break;
    case EventType::kMtoTEnd:
      return std::unique_ptr<Event>{new MtoTEnd{}};
      break;
    case EventType::kLtoTBegin:
      return std::unique_ptr<Event>{new LtoTBegin{}};
      break;
    case EventType::kLtoTEnd:
      return std::unique_ptr<Event>{new LtoTEnd{}};
      break;
    case EventType::kStoMBegin:
      return std::unique_ptr<Event>{new StoMBegin{}};
      break;
    case EventType::kStoMEnd:
      return std::unique_ptr<Event>{new StoMEnd{}};
      break;
    case EventType::kStoLBegin:
      return std::unique_ptr<Event>{new StoLBegin{}};
      break;
    case EventType::kStoLEnd:
      return std::unique_ptr<Event>{new StoLEnd{}};
      break;
    case EventType::kELCOBegin:
      return std::unique_ptr<Event>{new ELCOBegin{}};
      break;
    case EventType::kELCOEnd:
      return std::unique_ptr<Event>{new ELCOEnd{}};
      break;
    case EventType::kMtoMBegin:
      return std::unique_ptr<Event>{new MtoMBegin{}};
      break;
    case EventType::kMtoMEnd:
      return std::unique_ptr<Event>{new MtoMEnd{}};
      break;
    case EventType::kMtoLBegin:
      return std::unique_ptr<Event>{new MtoLBegin{}};
      break;
    case EventType::kMtoLEnd:
      return std::unique_ptr<Event>{new MtoLEnd{}};
      break;
    case EventType::kLtoLBegin:
      return std::unique_ptr<Event>{new LtoLBegin{}};
      break;
    case EventType::kLtoLEnd:
      return std::unique_ptr<Event>{new LtoLEnd{}};
      break;
    case EventType::kMtoIBegin:
      return std::unique_ptr<Event>{new MtoIBegin{}};
      break;
    case EventType::kMtoIEnd:
      return std::unique_ptr<Event>{new MtoIEnd{}};
      break;
    case EventType::kItoIBegin:
      return std::unique_ptr<Event>{new ItoIBegin{}};
      break;
    case EventType::kItoIEnd:
      return std::unique_ptr<Event>{new ItoIEnd{}};
      break;
    case EventType::kItoLBegin:
      return std::unique_ptr<Event>{new ItoLBegin{}};
      break;
    case EventType::kItoLEnd:
      return std::unique_ptr<Event>{new ItoLEnd{}};
      break;
    case EventType::kZeroRef:
      return std::unique_ptr<Event>{new ZeroRef{}};
      break;
    case EventType::kUnknown:
      return std::unique_ptr<Event>{nullptr};
      break;
  }

  return std::unique_ptr<Event>{nullptr};
}


const std::string &DASHMMEvent::event_class() const {
  return kEventClass;
}


const std::string &StoTBegin::event_type() const {
  return kStoTBeginType;
}

const std::string &StoTEnd::event_type() const {
  return kStoTEndType;
}

const std::string &MtoTBegin::event_type() const {
  return kMtoTBeginType;
}

const std::string &MtoTEnd::event_type() const {
  return kMtoTEndType;
}

const std::string &LtoTBegin::event_type() const {
  return kLtoTBeginType;
}

const std::string &LtoTEnd::event_type() const {
  return kLtoTEndType;
}

const std::string &StoMBegin::event_type() const {
  return kStoMBeginType;
}

const std::string &StoMEnd::event_type() const {
  return kStoMEndType;
}

const std::string &StoLBegin::event_type() const {
  return kStoLBeginType;
}

const std::string &StoLEnd::event_type() const {
  return kStoLEndType;
}

const std::string &ELCOBegin::event_type() const {
  return kELCOBeginType;
}

const std::string &ELCOEnd::event_type() const {
  return kELCOEndType;
}

const std::string &MtoMBegin::event_type() const {
  return kMtoMBeginType;
}

const std::string &MtoMEnd::event_type() const {
  return kMtoMEndType;
}

const std::string &MtoLBegin::event_type() const {
  return kMtoLBeginType;
}

const std::string &MtoLEnd::event_type() const {
  return kMtoLEndType;
}

const std::string &LtoLBegin::event_type() const {
  return kLtoLBeginType;
}

const std::string &LtoLEnd::event_type() const {
  return kLtoLEndType;
}

const std::string &MtoIBegin::event_type() const {
  return kMtoIBeginType;
}

const std::string &MtoIEnd::event_type() const {
  return kMtoIEndType;
}

const std::string &ItoIBegin::event_type() const {
  return kItoIBeginType;
}

const std::string &ItoIEnd::event_type() const {
  return kItoIEndType;
}

const std::string &ItoLBegin::event_type() const {
  return kItoLBeginType;
}

const std::string &ItoLEnd::event_type() const {
  return kItoLEndType;
}

const std::string &ZeroRef::event_type() const {
  return kZeroRefType;
}


std::unique_ptr<Event> StoTBegin::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new StoTBegin{rval}};
}

std::unique_ptr<Event> StoTEnd::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new StoTEnd{rval}};
}

std::unique_ptr<Event> MtoTBegin::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new MtoTBegin{rval}};
}

std::unique_ptr<Event> MtoTEnd::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new MtoTEnd{rval}};
}

std::unique_ptr<Event> LtoTBegin::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new LtoTBegin{rval}};
}

std::unique_ptr<Event> LtoTEnd::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new LtoTEnd{rval}};
}

std::unique_ptr<Event> StoMBegin::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new StoMBegin{rval}};
}

std::unique_ptr<Event> StoMEnd::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new StoMEnd{rval}};
}

std::unique_ptr<Event> StoLBegin::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new StoLBegin{rval}};
}

std::unique_ptr<Event> StoLEnd::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new StoLEnd{rval}};
}

std::unique_ptr<Event> ELCOBegin::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new ELCOBegin{rval}};
}

std::unique_ptr<Event> ELCOEnd::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new ELCOEnd{rval}};
}

std::unique_ptr<Event> MtoMBegin::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new MtoMBegin{rval}};
}

std::unique_ptr<Event> MtoMEnd::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new MtoMEnd{rval}};
}

std::unique_ptr<Event> MtoLBegin::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new MtoLBegin{rval}};
}

std::unique_ptr<Event> MtoLEnd::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new MtoLEnd{rval}};
}

std::unique_ptr<Event> LtoLBegin::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new LtoLBegin{rval}};
}

std::unique_ptr<Event> LtoLEnd::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new LtoLEnd{rval}};
}

std::unique_ptr<Event> MtoIBegin::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new MtoIBegin{rval}};
}

std::unique_ptr<Event> MtoIEnd::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new MtoIEnd{rval}};
}

std::unique_ptr<Event> ItoIBegin::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new ItoIBegin{rval}};
}

std::unique_ptr<Event> ItoIEnd::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new ItoIEnd{rval}};
}

std::unique_ptr<Event> ItoLBegin::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new ItoLBegin{rval}};
}

std::unique_ptr<Event> ItoLEnd::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new ItoLEnd{rval}};
}

std::unique_ptr<Event> ZeroRef::read_from_file(FILE *fd) const {
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
      throw std::runtime_error("error during read.");
    }
  }
  return std::unique_ptr<Event>{new ZeroRef{rval}};
}



} // traceutils::dashmm

} // traceutils
