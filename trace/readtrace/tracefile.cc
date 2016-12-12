#include "tracefile.h"

#include <exception>
#include <string>
#include <vector>

#include "allevents.h"
#include "traceevent.h"


namespace trace {

namespace {


  std::vector<std::string> split_once(const std::string &orig,
                                      const std::string &delim) {
    std::vector<std::string> retval{};

    auto floc = orig.find(delim, 0);
    if (floc == 0) {
      retval.push_back("");
      retval.push_back(orig.substr(floc + delim.size()));
    } else if (floc != std::string::npos) {
      retval.push_back(orig.substr(0, floc));
      retval.push_back(orig.substr(floc + delim.size()));
    } else {
      retval.push_back(orig);
      retval.push_back("");
    }

    return retval;
  }


  std::vector<std::string> split(const std::string &orig,
                                 const std::string &delim) {
    std::vector<std::string> retval{};

    std::string::size_type iter{0};
    auto delim_len = delim.size();

    while (iter != std::string::npos) {
      auto floc = orig.find(delim, iter);
      if (floc == iter) {
        iter = floc + delim_len;
      } else if (floc != std::string::npos) {
        auto piece_len = floc - iter;
        retval.push_back(orig.substr(iter, piece_len));
        iter = floc + delim_len;
      } else {
        retval.push_back(orig.substr(iter));
        iter = floc;
      }
    }

    return retval;
  }


  std::unique_ptr<Event> proto_from_class_and_name(EventClass evt_class,
                                                   const std::string &ename) {
    switch (evt_class) {
      case EventClass::kParcel:
        return std::unique_ptr<Event>{nullptr};
        break;
      case EventClass::kNetwork:
        return network::prototype_from_name(ename);
        break;
      case EventClass::kSched:
      case EventClass::kLCO:
      case EventClass::kProcess:
      case EventClass::kMemory:
      case EventClass::kTrace:
      case EventClass::kGAS:
      case EventClass::kCollective:
      case EventClass::kUnknown:
        return std::unique_ptr<Event>{nullptr};
        break;
    }

    return std::unique_ptr<Event>{nullptr};
  }


} // anonymous


File::File(const std::string &fname)
    : ifd_{nullptr}, proto_{nullptr}, worker_{-1}, locality_{-1} {
  ifd_ = fopen(fname.c_str(), "rb");
  if (ifd_ == nullptr) {
    // TODO: add some more meaning to this error - that is, copy out the
    // perror() message
    throw std::runtime_error("Unable to open file");
  }

  try {
    interpret_filename(fname);
  } catch (...) {
    // Some exception occurred, pass it along, after cleaning up resources
    fclose(ifd_);
    throw;
  }
}


File::~File() {
  fclose(ifd_);
}


// This will pass along the exception from the read
// This returns null when done
std::unique_ptr<Event> File::next() {
  return proto_->read_from_file(ifd_);
}


// This will possibly throw std::runtime_error
void File::interpret_filename(const std::string &fname) {
  auto pathparts = split(fname, "/");
  auto finalpathpart = pathparts[pathparts.size() - 1];

  auto parts = split(finalpathpart, ".");
  if (parts.size() != 6) {
    throw std::runtime_error("Filename does not match expectations");
  }

  try {
    locality_ = std::stoi(parts[0]);
  } catch (std::invalid_argument &invarg) {
    throw std::runtime_error("unable to convert locality");
  } catch (std::out_of_range &oorange) {
    throw std::runtime_error("unable to convert locality");
  }

  try {
    worker_ = std::stoi(parts[1]);
  } catch (std::invalid_argument &invarg) {
    throw std::runtime_error("unable to convert worker");
  } catch (std::out_of_range &oorange) {
    throw std::runtime_error("unable to convert worker");
  }

  auto c_and_e = split_once(parts[4], "_");
  if (c_and_e.size() != 2 || !(c_and_e[0].size() && c_and_e[1].size())) {
    throw std::runtime_error("error in parsing filename for event data");
  }

  EventClass evt_class = event_class_from_name(c_and_e[0]);
  proto_ = proto_from_class_and_name(evt_class, c_and_e[1]);
  if (proto_ == nullptr) {
    throw std::runtime_error("Unknown event class or type");
  }
}



} // trace