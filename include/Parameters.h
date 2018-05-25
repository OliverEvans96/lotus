#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "Utils.h"
#include <yaml.h>

using namespace std;

typedef map<string, vector<string> > StrVecMap;

///////////////////////////
// High level parameters //
///////////////////////////

struct Options
{
  StrVecMap yamlMap;
  char* configPath;

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
  string inLoc;
  string outLoc;

  vector<int> liquidTypes, solidTypes;

  StrVecMap parseYaml(char* configPath);
  bool mapHasKey(StrVecMap yamlMap, string key);
  void printYamlMap(StrVecMap yamlMap);
  void fromString(string optionString, bool &option);
  void fromString(string optionString, int &option);
  void fromString(string optionString, double &option);
  void fromString(string optionString, string &option);

  template <typename T>
  void unsafeParseOption(string optionName, T &option);
  template <typename T>
  void unsafeParseOption(string optionName, vector<T> &optionVec);
  template <typename T>
  void parseDefaultOption(string optionName, T &option, T defaultValue);
  template <typename T>
  void parseRequiredOption(string optionName, T &option);
  void readConfig(char* configPath);

  template <typename T>
  void printOption(string optionName, T option);
  template <typename T>
  void printOption(string optionName, vector<T> option);
  void printOptions();
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
  char *configPath;
  Options options;

  void parseArgs(int argc, char* argv[]);
  void print();
};

#endif
