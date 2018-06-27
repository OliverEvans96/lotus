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
  char configPath[256];

  string geometry; // spherical or cylindrical
  int stepsPerFrame;
  bool skipToEnd; // Jump to last frame
  bool trackMonoAtoms; // Follow a few atoms that end up in the monolayer
  bool saveImages; // Produce images or not
  bool saveROOT; // Produce ROOT figures or not
  bool plotHist; // Whether to draw 2D cylindrical histogram
  bool plotDipole;
  bool plotVr;
  bool plotDensity;
  bool plotAllTogether;
  bool verbose; // Enable debuguging print statments
  bool onlyFindInterface; // Locate monolayer in first frame and exit
  string dumpfile;
  string datafile;
  string outLoc;
  int waterBondType;
  double dz;
  double dv;
  double dens_min;
  double dens_max;
  double plot_rmax;
  double plot_zmax;
  double plot_aspect;
  int plot_width;
  double expectedLiquidDensity;

  vector<int> liquidTypes, solidTypes;

  StrVecMap parseYaml(const char* configPath);
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
  void readConfig(const char* configPath);

  template <typename T>
  void printOption(string optionName, T option);
  template <typename T>
  void printOption(string optionName, vector<T> option);
  void printOptions();
};

struct CommandLineParser
{
  CommandLineParser(int argc, const char* argv[]);
  char configPath[256];
  Options options;

  void parseArgs(int argc, const char* argv[]);
  void print();
};

#endif
