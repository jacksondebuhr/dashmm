#include <cstdio>

#include <algorithm>
#include <exception>
#include <memory>
#include <string>
#include <vector>

#include "locality.h"
#include "plotter.h"
#include "segmentids.h"
#include "traceevent.h"
#include "tracefile.h"
#include "trace.h"
#include "utilization.h"
#include "worker.h"


void file_report(const char *fname, const traceutils::File &file) {
  fprintf(stdout, "File '%s' opened...\n", fname);
  fprintf(stdout, "  Locality: %d\n", file.locality());
  fprintf(stdout, "  Worker: %d\n", file.worker());
  fprintf(stdout, "  Event Class: '%s'\n", file.event_class().c_str());
  fprintf(stdout, "  Event Type: '%s'\n", file.event_type().c_str());
}


void do_output(FILE *ofd, traceutils::iter_t begin, traceutils::iter_t end,
               uint64_t span_start) {
  for (auto i = begin; i != end; ++i) {
    fprintf(ofd, "%.12lg %s\n", ((*i)->stamp() - span_start) / 1.0e6,
            (*i)->event_type().c_str());
  }
}

void output_window(const traceutils::window_t &window, uint64_t span_start) {
  FILE *ofd = fopen("eventlog.txt", "w");

  for (auto i = window.begin(); i != window.end(); ++i) {
    for (auto j = i->second.begin(); j != i->second.end(); ++j) {
      do_output(ofd, j->second.first, j->second.second, span_start);
    }
  }

  fclose(ofd);
}


int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stdout, "Usage: %s <trace file name>\n", argv[0]);
    return 0;
  }

  traceutils::Trace runtrace;

  try {
    for (int which = 1; which < argc; ++which) {
      traceutils::File in_file{std::string(argv[which])};
      file_report(argv[which], in_file);
      runtrace.add_file(in_file);
    }
    runtrace.find_and_reset_zero();
    runtrace.finalize();

    uint64_t minns = runtrace.min_ns();
    uint64_t maxns = runtrace.max_ns();
    fprintf(stdout, "\nFiles contained %lu events\n", runtrace.num_events());
    fprintf(stdout, "Time window: [%lg %lg] (ms)\n",
            minns / 1.0e6, maxns / 1.0e6);

    // Get a window of events from 1s to 2s
    //auto window = runtrace.window(minns, maxns);
    //output_window(window, runtrace.min_ns());

    //auto coverage = coverage_of_segment_type(window, 1,
    //                                         minns, maxns);

    //for (auto i = coverage.begin(); i != coverage.end(); ++i) {
    //  for (auto j = i->second.begin(); j != i->second.end(); ++j) {
    //    fprintf(stdout, "coverage of (%d,%d): %lg\n", i->first, j->first,
    //            (coverage[i->first])[j->first]);
    //  }
    //}

    /*
    for (int segment = traceutils::segment::kNetworkProgress;
         segment <= traceutils::segment::kDASHMMItoL; ++segment) {
      traceutils::Plotter maker{runtrace};
      auto dt = maxns - minns;
      std::string segtype = traceutils::segment::name(segment);

      maker(1000, 512, 0, minns, maxns, segtype + "full.png", segment, 0);
      maker(1000, 512, 0, minns, minns + dt * 0.25, segtype + "first.png",
            segment, 0);
      maker(1000, 512, 0, minns + dt * 0.25, minns + dt * 0.5,
            segtype + "second.png", segment, 0);
      maker(1000, 512, 0, minns + dt * 0.5, minns + dt * 0.75,
            segtype + "third.png", segment, 0);
      maker(1000, 512, 0, minns + dt * 0.75, maxns,
            segtype + "fourth.png", segment, 0);
    }
    */

    std::vector<int> segs{};
    for (int seg = traceutils::segment::kDASHMMStoT;
         seg <= traceutils::segment::kDASHMMItoL; ++seg) {
      segs.push_back(seg);
    }

    traceutils::Utilization util{runtrace};
    util("test.txt", minns, maxns, 100, segs);


  } catch (std::runtime_error &err) {
    fprintf(stderr, "Exception: %s\n", err.what());
  } catch (std::invalid_argument &err) {
    fprintf(stderr, "Exception: %s\n", err.what());
  } catch (...) {
    fprintf(stderr, "Who knows...\n");
  }


  // Testing out the thing somewhat
  //pngwriter testimage{100, 100, 1.0, "checker.png"};
  //testimage.setcompressionlevel(0);
  //testimage.filledsquare(10, 10, 50, 50, 1.0, 0.0, 0.0);
  //testimage.close();

  return 0;
}
