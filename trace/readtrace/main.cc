#include <cstdio>

#include <algorithm>
#include <exception>
#include <memory>
#include <string>
#include <vector>

#include "locality.h"
#include "plotter.h"
#include "segmentids.h"
#include "sqlite3out.h"
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

    /*
    auto window = runtrace.window(minns, maxns);
    std::vector<int> segs{};
    for (int seg = traceutils::segment::kDASHMMStoT;
         seg <= traceutils::segment::kDASHMMItoL; ++seg) {
      segs.push_back(seg);
    }

    traceutils::Utilization util{runtrace};
    util("test.txt", minns, maxns, 100, segs);
    //*/

    // now go ahead and put it out as an SQLite database
    //*
    traceutils::SQLiteWriter out{"trace.db", runtrace};
    out.write();
    //*/


  } catch (std::runtime_error &err) {
    fprintf(stderr, "Exception: %s\n", err.what());
  } catch (std::invalid_argument &err) {
    fprintf(stderr, "Exception: %s\n", err.what());
  } catch (...) {
    fprintf(stderr, "Who knows...\n");
  }

  return 0;
}
