#include "Atoms.h"
#include "Readers.h"
#include "catch.hpp"

using namespace std;

TEST_CASE("Readers", "[lotus]") {
  AtomArray *atomArrayPtr;
  SimData *simDataPtr;
  int argc = 2;
  char* argv[2] = {"test", "../test/data/test_config.yaml"};
  CommandLineParser commandLineParser(argc, argv);
  Options options = commandLineParser.options;
  FrameReader frameReader(options, atomArrayPtr, simDataPtr);
}
