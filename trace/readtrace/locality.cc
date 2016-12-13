#include "locality.h"

#include <exception>


namespace traceutils {


Locality::Locality(File &stream) : locality_{stream.locality()}, workers_{} {
  Worker to_add{stream};
  workers_[stream.worker()] = std::move(to_add);
}


size_t Locality::num_events() const {
  size_t retval{};
  for (auto i = workers_.begin(); i != workers_.end(); ++i) {
    retval += i->second.num_events();
  }
  return retval;
}


void Locality::add_file(File &stream) {
  if (stream.locality() != locality_) {
    throw std::invalid_argument("File does not match locality");
  }
  auto wkr = workers_.find(stream.worker());
  if (wkr != workers_.end()) {
    wkr->second.add_file(stream);
  } else {
    Worker to_add{stream};
    workers_[stream.worker()] = std::move(to_add);
  }
}


} // traceutils
