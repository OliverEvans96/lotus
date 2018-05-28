#ifndef READERS_H
#define READERS_H

#include <cstdio>
#include <map>
#include <cctype>

#include "Utils.h"
#include "MDBase.h"
#include "Atoms.h"
#include "Droplet.h"
#include "Time.h"

using namespace std;

//////////////////
// Input Stream //
//////////////////

struct InputStream {
  int lineNum;
  string filename;
  ifstream stream;

  InputStream();
  InputStream(string _filename);
  ~InputStream();
  void open(string _filename);
  void verifyStream();
  void skipLines(int numLines);
  void skipWhitespace();
  bool search(string term);
  string search(vector<string> terms);
  bool searchLine(string term);
  string searchLine(vector<string> terms);
  bool nextLineBlank();
  string peekLine();
};

//Split a string into a string vector of words
vector<string> strSplit(string str);

//Get index and coordinates from string array of words
void strToData(double *coords,double *velocities,double& dipole,string line);


/////////////
// Readers //
/////////////

// These read from an InputStream, but don't actually open the files.
class HeaderReader {
  Options options;
  InputStream* inputStreamPtr;
  Timestep* timestepPtr;

 public:
  int* lineNumPtr;
  void setContext(Options options, InputStream* _inputStreamPtr, Timestep* _timestepPtr, int* _lineNumPtr);
  void readHeader();
};

class LineReader {
  Options options;
  string line;
  string input;
  int* atomNumPtr;
  int* lineNumPtr;
  InputStream* inputStreamPtr;
  vector<string> strSplit(string str);
  void strToData(double* coords, double* velocities, double& dipole, string line);

 public:
  Atom atom;
  void setContext(Options options, InputStream *_inputStreamPtr, int* _atomNumPtr, int* _lineNumPtr);
  void readLine();
};

class TimestepReader {
  Options options;
  LineReader lineReader;
  HeaderReader headerReader;
  Timestep* timestepPtr;
  AtomArray* atomArrayPtr;
  InputStream* inputStreamPtr;
  SimData* simDataPtr;

 public:
  int atomNum;
  int lineNum;

  TimestepReader();
  void setContext(Options _options, InputStream* _inputStreamPtr, AtomArray* _atomArrayPtr, Timestep* _timestepPtr, SimData* _simDataPtr);
  void resetAtomCounter();
  void readTimestep();
};

class FrameReader {
 public:

  Options options;
  TimestepReader timestepReader;
  AtomArray* atomArrayPtr;
  Timestep timestep;
  SimData* simDataPtr;
  InputStream inputStream;
  Frame frame;
  int stepsPerFrame;

  FrameReader();
  FrameReader(Options options, AtomArray* _atomArrayPtr, SimData* _simDataPtr);
  void setContext(Options options, AtomArray* _atomArrayPtr, SimData* _simDataPtr);
  void openStream();
  void updateFrame();
  void readFrame();

  void countAtoms();
  void countSteps();
};

class InitialTimestepReader {
  Options options;
 public:
  InitialTimestepReader(Options options, AtomArray* _atomArrayPtr, SimData* _simDataPtr);
  void setContext(Options options, AtomArray* _atomArrayPtr, SimData* _simDataPtr);
  string initLoc;
  TimestepReader timestepReader;
  AtomArray* atomArrayPtr;
  Timestep* timestepPtr;
  InputStream inputStream;
  Timestep emptyTimestep;
  SimData* simDataPtr;
  Atom atom;

  bool fileExists();
  void openStream(Options _options);
  void readFromFile();
  void readFromStream();
  void writeToFile();
  void readInitialTimestep();
};

///////////////////////
// Top-level Readers //
///////////////////////

// Read data file for all Masses & liquid Bonds
class DatafileReader {
  Options options;
  InputStream inputStream;
  vector<int> liquidTypes;
  int HType, OType;
  ifstream *streamPtr;
  SimData *simDataPtr;
  int numAtoms;

  void readNumAtoms();
  void readMasses();
  void readWaterBonds();

 public:
  DatafileReader(Options _options, SimData &simData);
  void openStream();
  void read();
};

struct DumpfileReader {
  Options options;
  FrameReader frameReader;
  SimData *simDataPtr;
  AtomArray *atomArrayPtr;

  DumpfileReader(Options _options, SimData &simData, AtomArray &atomArray);
  void countAtoms();
  void countSteps();
  void readFrame();
  bool good();
};


#endif
