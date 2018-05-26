#include "Parameters.h"
#include "MDBase.h"
#include "Atoms.h"
#include "Readers.h"
#include "catch.hpp"

using namespace std;

const int NUM_ATOMS = 513699;

TEST_CASE("Readers", "[lotus]") {
  int argc = 2;
  char* argv[2] = {"test", "../test/data/test_config.yaml"};
  CommandLineParser commandLineParser(argc, argv);
  Options options = commandLineParser.options;
  SimData simData(options);
  AtomArray atomArray(simData);

  REQUIRE(simData.numAtoms == NUM_ATOMS);

  FrameReader frameReader(options, &atomArray, &simData);
  REQUIRE(atomArray.numAtoms == NUM_ATOMS);

  //frameReader.readFrame();
}
