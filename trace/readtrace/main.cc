#include <cstdio>

#include <algorithm>
#include <exception>
#include <memory>
#include <string>
#include <vector>

#include "locality.h"
#include "traceevent.h"
#include "tracefile.h"
#include "trace.h"
#include "worker.h"


void file_report(const char *fname, const traceutils::File &file) {
  fprintf(stdout, "File '%s' opened...\n", fname);
  fprintf(stdout, "  Locality: %d\n", file.locality());
  fprintf(stdout, "  Worker: %d\n", file.worker());
  fprintf(stdout, "  Event Class: '%s'\n", file.event_class().c_str());
  fprintf(stdout, "  Event Type: '%s'\n", file.event_type().c_str());
}


void do_output(const std::vector<std::unique_ptr<traceutils::Event>> &events) {
  FILE *ofd = fopen("eventlog.txt","w");
  for (size_t i = 0; i < events.size(); ++i) {
    fprintf(ofd, "%lu %s\n", events[i]->stamp(),
            events[i]->event_type().c_str());
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

    fprintf(stdout, "\nFiles contained %lu events\n", runtrace.num_events());
    fprintf(stdout, "Time window: [%lg %lg] (ms)\n",
            runtrace.min_ns() / 1.0e6, runtrace.max_ns() / 1.0e6);

    fprintf(stdout, "Finalizing..."); fflush(stdout);
    runtrace.finalize();
    fprintf(stdout, "done\n"); fflush(stdout);

    //do_output(events);

    // Get a window of events from 1s to 2s
    auto window = runtrace.window(1000000000, 2000000000);

  } catch (std::runtime_error &err) {
    fprintf(stderr, "Exception: %s\n", err.what());
  } catch (...) {
    fprintf(stderr, "Some other exception...\n");
  }

  return 0;
}
