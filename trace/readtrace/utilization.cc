#include "utilization.h"

#include <cassert>
#include <cstdio>


namespace traceutils {


namespace {

  double sum_of_cover(const cover_t &c) {
    double retval{0.0};
    int count{0};

    for (auto i = c.begin(); i != c.end(); ++i) {
      for (auto j = i->second.begin(); j != i->second.end(); ++j) {
        ++count;
        retval += j->second;
      }
    }

    assert(count);

    return retval / count;
  }

} // {anonymous}


void Utilization::operator()(const std::string &fname, uint64_t t0, uint64_t t1,
                             int samples, const std::vector<int> segments) {
  // The first index is time, the second index is column
  std::vector<std::vector<double>> output(samples);
  for (int i = 0; i < samples; ++i) {
    // +2 because we have the sum and the timestamp
    output[i].reserve(segments.size() + 2);
  }

  // Compute some one-time things
  double dx = (double)(t1 - t0) / samples;

  for (int s = 0; s < samples; ++s) {
    uint64_t tend = (uint64_t)(dx * (s + 1) + 0.5) + t0;
    uint64_t tstart = (uint64_t)(dx * s + 0.5) + t0;
    output[s][0] = dx * (s + 0.5);

    auto window = trace_.window(tstart, tend);
    for (size_t i = 0; i < segments.size(); ++i) {
      auto coverage = coverage_of_segment_type(window, segments[i],
                                               tstart, tend);
      // +1 to skip the timestamp
      output[s][i + 1] = sum_of_cover(coverage);
    }
  }

  for (int s = 0; s < samples; ++s) {
    double total{0.0};
    for (size_t i = 0; i < segments.size(); ++i) {
      total += output[s][i + 1];
    }
    output[s][segments.size() + 1] = total;
  }

  // Now open the file and print it out
  FILE *ofd = fopen(fname.c_str(), "w");
  assert(ofd);

  for (int i = 0; i < samples; ++i) {
    for (size_t j = 0; j < segments.size() + 2; ++j) {
      fprintf(ofd, "%lg ", output[i][j]);
    }
    fprintf(ofd, "\n");
  }

  fclose(ofd);
}


} // traceutils