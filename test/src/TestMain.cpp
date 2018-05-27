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
  REQUIRE(options.configPath == argv[1]);
  REQUIRE(options.liquidTypes[0] == 4);
  REQUIRE(options.liquidTypes[1] == 5);
  REQUIRE(options.solidTypes[0] == 1);
  REQUIRE(options.solidTypes[1] == 2);
  REQUIRE(options.solidTypes[2] == 3);
  REQUIRE(options.geometry == "spherical");
  REQUIRE(options.stepsPerFrame == 3);
  REQUIRE(options.dumpfile == "../test/data/20A_atom2_first8");
  REQUIRE(options.datafile == "../test/data/lammps_noZperiod_3A.dat");
  REQUIRE(options.outLoc == "../test/run");
  REQUIRE(options.skipToEnd == false);
  REQUIRE(options.trackMonoAtoms == false);
  REQUIRE(options.saveImages == false);
  REQUIRE(options.plotHist == false);
  REQUIRE(options.plotDipole == false);
  REQUIRE(options.plotVr == false);
  REQUIRE(options.plotDensity == false);
  REQUIRE(options.plotAllTogether == false);
  REQUIRE(options.debugOutput == false);
  REQUIRE(options.onlyFindInterface == false);
  SimData simData(options);
  Droplet droplet(options);
  DatafileReader datafileReader(options, simData);
  REQUIRE(simData.masses.size() == 5);
  REQUIRE(simData.masses[1] == 26.981540);
  REQUIRE(simData.masses[2] == 15.994915);
  REQUIRE(simData.masses[3] == 1.007825);
  REQUIRE(simData.masses[4] == 15.9994);
  REQUIRE(simData.masses[5] == 1.00794);

  REQUIRE(simData.numAtoms == NUM_ATOMS);
  REQUIRE(simData.waterBonds.size() == NUM_WATER);

  DumpfileReader dumpfileReader(options, simData, droplet);
  REQUIRE(simData.numSteps == NUM_STEPS);

  // Time loop
  while(dumpfileReader.good()) {
    dumpfileReader.readFrame();
  }
}

