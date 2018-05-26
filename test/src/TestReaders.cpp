#include "Parameters.h"
#include "MDBase.h"
#include "Atoms.h"
#include "Readers.h"
#include "catch.hpp"

using namespace std;

const int NUM_ATOMS = 513699;
const int NUM_STEPS = 8;

TEST_CASE("Readers", "[lotus]") {
  int argc = 2;
  char* argv[2] = {"test", "../test/data/test_config.yaml"};
  CommandLineParser commandLineParser(argc, argv);
  Options options = commandLineParser.options;
  SimData simData(options);

  DatafileReader datafileReader(options, simData);
  REQUIRE(simData.numAtoms == NUM_ATOMS);

  DumpfileReader dumpfileReader(options, simData);
  REQUIRE(simData.numSteps == NUM_STEPS);

  //frameReader.readFrame();
}
