#ifndef __TRACEUTILS_LOCALITY_H__
#define __TRACEUTILS_LOCALITY_H__


#include <map>

#include "tracefile.h"
#include "worker.h"


namespace traceutils {


class Locality {
 public:
  Locality(File &stream);
  Locality(int loc = -1) : locality_{loc}, workers_{} { }
  Locality(Locality &&other) = default;

  // Make sure default move operator exists
  Locality &operator=(Locality &&other) = default;


  int locality() const {return locality_;}
  // The number of workers currently.
  int num_workers() const {return workers_.size();}
  // The maximum id of any workers.
  int max_worker() const {return workers_.rbegin()->first;}

  size_t num_events() const;

  // This could throw std::runtime_error or std::invalid_argument
  void add_file(File &stream);

  // TODO: more complicated queries - not really sure what I need yet
 private:
  int locality_;
  std::map<int, Worker> workers_;
};


} // traceutils


#endif // __TRACEUTILS_LOCALITY_H__