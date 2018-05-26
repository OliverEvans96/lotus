#include "Parameters.h"
#include "catch.hpp"

using namespace std;

// int main(int argc, char* argv[]) {
//   CommandLineParser commandLineParser(argc, argv);
//   commandLineParser.print();
//   return 0;
// }

TEST_CASE("Parameters", "[lotus]") {
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
  REQUIRE(options.stepsPerFrame == 3);
  REQUIRE(options.inLoc == "../test/data/20A_atom2_first8");
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
}
