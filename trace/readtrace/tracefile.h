#ifndef __TRACE_FILE_H__
#define __TRACE_FILE_H__


#include <cstdio>

#include <memory>
#include <string>

#include "traceevent.h"


namespace trace {


class File {
 public:
  // Should we have this emit exceptions for bad formats and things
  // Probably. It will be interesting to do exceptions all the way once to
  // see how that all hangs together.
  File(const std::string &fname);
  ~File();

  // basically I want the interface to be as a factory for events
  // and the event type will depend on the file and so on.
  //
  // This is going to emit an exception to indicate errors. It is just going to
  // be std::runtime_error with a relevant error message.
  //
  // This will return a null pointer if the eof is found at the start of a
  // new record
  std::unique_ptr<Event> next();

  // These are all detected from the filename
  std::string event_class() const {return proto_->event_class();}
  std::string event_type() const {return proto_->event_type();}
  int worker() const noexcept {return worker_;}
  int locality() const noexcept {return locality_;}

  // Allow something to close it early
  void close();

 private:
  // utility function to pull name apart to get class and event of the file
  // as well as the worker and the locality
  void interpret_filename(const std::string &fname);

  // utility function to skip the header information
  void check_and_skip_header();

  FILE *ifd_;
  std::unique_ptr<Event> proto_;
  int worker_;
  int locality_;
};


} // trace


#endif // __TRACE_FILE_H__