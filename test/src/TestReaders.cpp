#include "Parameters.h"
#include "MDBase.h"
#include "Atoms.h"
#include "Readers.h"
#include "catch.hpp"

using namespace std;

const int NUM_ATOMS = 513699;
const int NUM_STEPS = 8;
const int NUM_WATER = 1105;

TEST_CASE("Readers", "[lotus]") {
  int argc = 2;
  char* argv[2] = {"test", "../test/data/test_config.yaml"};
  CommandLineParser commandLineParser(argc, argv);
  Options options = commandLineParser.options;
  SimData simData(options);

  DatafileReader datafileReader(options, simData);
  REQUIRE(simData.masses.size() == 5);
  REQUIRE(simData.masses[1] == 26.981540);
  REQUIRE(simData.masses[2] == 15.994915);
  REQUIRE(simData.masses[3] == 1.007825);
  REQUIRE(simData.masses[4] == 15.9994);
  REQUIRE(simData.masses[5] == 1.00794);

  REQUIRE(simData.numAtoms == NUM_ATOMS);
  REQUIRE(simData.waterBonds.size() == NUM_WATER);

  DumpfileReader dumpfileReader(options, simData);
  REQUIRE(simData.numSteps == NUM_STEPS);

  //dumpfileReader.readFrame();
}

