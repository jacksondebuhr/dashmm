#ifndef __TRACEUTILS_COUNTER_H__
#define __TRACEUTILS_COUNTER_H__


namespace traceutils {


class Counter {
 public:
  Counter() : c_{0} { }
  uint64_t next() {return c_++;}

  // Not really sure if this is worth having
  uint64_t curr() {return c_;}

 private:
  uint64_t c_;
};


}


#endif // __TRACEUTILS_COUNTER_H__