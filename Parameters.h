#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "Utils.h"

using namespace std;

///////////////////////////
// High level parameters //
///////////////////////////

struct Options
{
  bool skipToEnd; // Jump to last frame
  bool trackMonoAtoms; // Follow a few atoms that end up in the monolayer
  bool saveImages; // Produce images or not
  bool plotHist; // Whether to draw 2D cylindrical histogram
  bool plotDipole;
  bool plotVr;
  bool plotDensity;
  bool plotAllTogether;
  bool debugOutput; // Enable debuguging print statments
  bool onlyFindInterface; // Locate monolayer in first frame and exit

  void readOptions(string optionsFile);
  void print();
};

struct AnalysisParameters
{
  double rDensCyl; // Radius of cylinder used to calculate water density (z)
  double rBulkMax; // Maximum cylinder to use for radial binning
  double binVolume; // Constant volume for all cylindrical bins
};

struct CommandLineParser
{
  CommandLineParser(int argc, char* argv[]);
  string inLoc, outLoc, optionsFile;
  Options options;

  void parseArgs(int argc, char* argv[]);
  void print();
};

#endif
