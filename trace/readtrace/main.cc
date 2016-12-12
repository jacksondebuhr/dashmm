#include <cstdio>

#include <exception>
#include <memory>
#include <string>
#include <vector>

#include "traceevent.h"
#include "tracefile.h"


int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stdout, "Usage: %s <trace file name>\n", argv[0]);
    return 0;
  }

  // now try to create a TraceFile object
  try {
    trace::File in_file{std::string(argv[1])};

    fprintf(stdout, "File '%s' opened...\n", argv[1]);
    fprintf(stdout, "  Locality: %d\n", in_file.locality());
    fprintf(stdout, "  Worker: %d\n", in_file.worker());
    fprintf(stdout, "  Event Class: '%s'\n", in_file.event_class().c_str());
    fprintf(stdout, "  Event Type: '%s'\n", in_file.event_type().c_str());

    std::vector<std::unique_ptr<trace::Event>> events{};
    do {
      events.push_back(in_file.next());
    } while (events[events.size() - 1] != nullptr);
    events.pop_back();

    fprintf(stdout, "  File contained %lu events\n", events.size());

  } catch (std::runtime_error &err) {
    fprintf(stderr, "Exception: %s\n", err.what());
  } catch (...) {
    fprintf(stderr, "Some other exception...\n");
  }

  return 0;
}
