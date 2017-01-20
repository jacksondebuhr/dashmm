#include "tracefile.h"

#include <cstring>

#include <exception>
#include <vector>

#include "allevents.h"


namespace traceutils {

namespace {

  char file_magic_number[] = "HPXnpy";

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
        return dashmm::prototype_from_name(ename);
        break;
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
    check_and_skip_header();
  } catch (...) {
    // Some exception occurred, pass it along, after cleaning up resources
    fclose(ifd_);
    throw;
  }
}


File::~File() {
  close();
}


// This will pass along the exception from the read
// This returns null when done
std::unique_ptr<Event> File::next() {
  return proto_->read_from_file(ifd_);
}


void File::close() {
  if (ifd_) {
    fclose(ifd_);
    ifd_ = nullptr;
  }
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


void File::check_and_skip_header() {
  // read in the first 6 bytes
  char magic[6];
  if (6 != fread(magic, sizeof(char), 6, ifd_)) {
    throw std::runtime_error("Error reading header information.");
  }
  if (strncmp(magic, file_magic_number, 6)) {
    throw std::runtime_error("Format error: magic number is wrong");
  }

  if (fseek(ifd_, 2, SEEK_CUR)) {
    throw std::runtime_error("Error reading header information.");
  }
  uint16_t header_size{0};
  if (1 != fread(&header_size, sizeof(header_size), 1, ifd_)) {
    throw std::runtime_error("Error reading header information.");
  }
  if (fseek(ifd_, header_size, SEEK_CUR)) {
    throw std::runtime_error("Error reading header information.");
  }
}


} // traceutils