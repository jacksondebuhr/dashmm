#ifndef __TRACEUTILS_PLOTTER_H__
#define __TRACEUTILS_PLOTTER_H__


#include <string>

#include <pngwriter.h>

#include "trace.h"


namespace traceutils {


class Plotter {
 public:
  Plotter(const Trace &trace) : trace_{trace} { }

  // note py is a guideline - it will be set to its real value by the
  // constructor.
  // TODO: eventially we add segment types and colors to this. For now, we
  //  are fixed with red and network progress
  void operator()(int px, int py, int padding, uint64_t t0, uint64_t t1,
                  const std::string &fname, int segment, int compress = -1);

 private:
  struct image_data_t {
    image_data_t(int x, int y, int p, int h, uint64_t ti, uint64_t tf,
                 int nwork, const std::string &fname, int comp);

    int px;
    int py;
    int padding;
    int bar_height;
    int n_bars;
    uint64_t t0;
    uint64_t t1;
    pngwriter outfile;

    void set_bar_height();
    int bar_bottom(int idx);
    int bar_top(int idx);
    // TODO routine to get a certain bar's location
    void plot(const Trace &trace, int segment);
  };

  const Trace &trace_;
};


} // traceutils


#endif // __TRACEUTILS_PLOTTER_H__