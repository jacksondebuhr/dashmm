#include "plotter.h"


namespace traceutils {


void Plotter::operator()(int px, int py, int padding, uint64_t t0, uint64_t t1,
                         const std::string &fname, int segment, int compress) {
  image_data_t img{px, py, padding, 0, t0, t1, trace_.total_workers(),
                   fname, compress};
  img.plot(trace_, segment);
}



Plotter::image_data_t::image_data_t(
    int x, int y, int p, int h, uint64_t ti, uint64_t tf,
    int nwork, const std::string &fname, int comp) {
  px = x;
  py = y;
  padding = p;
  n_bars = nwork;
  set_bar_height();
  t0 = ti;
  t1 = tf;
  outfile = pngwriter(px, py, 1.0, fname.c_str());
  outfile.setcompressionlevel(comp);
}


void Plotter::image_data_t::set_bar_height() {
  double top = py - padding * (n_bars - 1);
  double bottom = n_bars;
  double x = top / bottom;
  bar_height = (int)(x + 0.5);
  py = bar_height * n_bars + padding * (n_bars - 1);
}


int Plotter::image_data_t::bar_bottom(int idx) {
  return (idx * (padding + bar_height));
}


int Plotter::image_data_t::bar_top(int idx) {
  return bar_bottom(idx) + bar_height;
}


void Plotter::image_data_t::plot(const Trace &trace, int segment) {
  // compute some one time things
  double dx = (double)(t1 - t0) / px;
  std::map<int, std::map<int, int>> wmap = trace.worker_map();
  int index = n_bars;
  for (auto i = wmap.begin(); i != wmap.end(); ++i) {
    for (auto j = i->second.begin(); j != i->second.end(); ++j) {
      j->second = --index;
    }
  }

  // loop over horizontal pixels
  for (int ix = 0; ix < px; ++ix) {
    uint64_t tstart = (uint64_t)(dx * ix + 0.5) + t0;
    uint64_t tend = (uint64_t)(dx * (ix + 1) + 0.5) + t0;
    auto window = trace.window(tstart, tend);
    auto coverage = coverage_of_segment_type(window, segment, tstart, tend);

    for (auto i = coverage.begin(); i != coverage.end(); ++i) {
      for (auto j = i->second.begin(); j != i->second.end(); ++j) {
        int idx = (wmap[i->first])[j->first];
        outfile.filledsquare(ix+1, bar_bottom(idx)+1, ix+1, bar_top(idx),
                             1.0, 1.0 - j->second, 1.0 - j->second);
      }
    }
  }

  outfile.close();
}


} // traceutils
