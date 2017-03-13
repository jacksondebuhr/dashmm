#include <cstdio>

#include "combiner.h"


// Print out a short usage statment
void print_usage(char *progname) {
  fprintf(stdout, "Usage: %s <file in 1> <file in 2> <file in 3> <file out>\n",
          progname);
}

int main(int argc, char **argv) {
  if (argc < 5) {
    print_usage(argv[0]);
    return -1;
  }

  try {
    Combiner fileone{argv[1]};

    { // We scope these to remove their resources as soon as possible
      Combiner filetwo{argv[2]};
      Combiner filethree{argv[3]};

      fileone.Combine(filetwo);
      filetwo.Combine(filethree);
    }

    fileone.Write(argv[4]);

  } catch (const std::ios_base::failure &e) {
    fprintf(stderr, "Bad news:\n%s\n", e.what());
  } catch (const std::runtime_error &e) {
    fprintf(stderr, "Watch out!\n%s\n", e.what());
  }

  fprintf(stdout,
          "Files: '%s', '%s' and '%s' successfully combined into '%s'\n",
          argv[1], argv[2], argv[3], argv[4]);

  return 0;
}