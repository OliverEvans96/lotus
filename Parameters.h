#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "Utils.h"

using namespace std;

///////////////////////////
// High level parameters //
///////////////////////////

struct CommandLineParser
{
  string inLoc, outLoc;
  void parse();
};

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
};

struct AnalysisParameters
{
  double rDensCyl; // Maximum cylinder to use for radial binning
  double rBulkMax; // ??
  double binVolume; // Constant volume for all cylindrical bins
};


#endif
