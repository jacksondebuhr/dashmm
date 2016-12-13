#include <cstdio>

#include <algorithm>
#include <exception>
#include <memory>
#include <string>
#include <vector>

#include "locality.h"
#include "traceevent.h"
#include "tracefile.h"
#include "worker.h"


void file_report(const char *fname, const trace::File &file) {
  fprintf(stdout, "File '%s' opened...\n", fname);
  fprintf(stdout, "  Locality: %d\n", file.locality());
  fprintf(stdout, "  Worker: %d\n", file.worker());
  fprintf(stdout, "  Event Class: '%s'\n", file.event_class().c_str());
  fprintf(stdout, "  Event Type: '%s'\n", file.event_type().c_str());
}


void do_output(const std::vector<std::unique_ptr<trace::Event>> &events) {
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

  try {
    trace::File first_file{std::string(argv[1])};
    file_report(argv[1], first_file);
    trace::Locality locality{first_file};

    for (int which = 2; which < argc; ++which) {
      trace::File in_file{std::string(argv[which])};
      file_report(argv[which], in_file);
      locality.add_file(in_file);
    }

    fprintf(stdout, "\n Files contained %lu events\n", locality.num_events());

    //do_output(events);

  } catch (std::runtime_error &err) {
    fprintf(stderr, "Exception: %s\n", err.what());
  } catch (...) {
    fprintf(stderr, "Some other exception...\n");
  }

  return 0;
}
