#ifndef READERS_H
#define READERS_H

#include <cstdio>

#include "Utils.h"
#include "MDBase.h"
#include "Atoms.h"
#include "Time.h"

using namespace std;

/////////////
// Readers //
/////////////

// These read from an InputStream, but don't actually open the files.


class HeaderReader {
  InputStream* inputStreamPtr;
  Timestep* timestepPtr;
  int* lineNumPtr;

 public:
  void setContext(InputStream* _inputStreamPtr, Timestep* _timestepPtr, int* lineNumPtr);
  void readHeader();
};

class LineReader {
  string line;
  string input;
  int* atomNumPtr;
  int* lineNumPtr;
  InputStream* inputStreamPtr;
  vector<string> strSplit(string str);
  void strToData(double* coords, double* velocities, double& dipole, string line);

 public:
  Atom atom;
  void setContext(InputStream *_inputStreamPtr, int* _atomNumPtr, int* _lineNumPtr);
  void readLine();
};

class TimestepReader {
  LineReader lineReader;
  HeaderReader headerReader;
  Timestep* timestepPtr;
  AtomArray* atomArrayPtr;
  InputStream* inputStreamPtr;
  SimData* simDataPtr;

 public:
  int atomNum;
  int lineNum;

  void setContext(InputStream* _inputStreamPtr, AtomArray* _atomArrayPtr, Timestep* _timestepPtr, SimData* _simDataPtr);
  void resetAtomCounter();
  void readTimestep();
};

class FrameReader {
 public:
  FrameReader(CommandLineParser commandLineParser, AtomArray* _atomArrayPtr, SimData* _simDataPtr);
  void setContext(CommandLineParser commandLineParser, AtomArray* _atomArrayPtr, SimData* _simDataPtr);


  TimestepReader timestepReader;
  AtomArray* atomArrayPtr;
  Timestep* timestepPtr;
  SimData* simDataPtr;
  InputStream inputStream;
  Frame frame;

  int stepsPerFrame;
  void openStream(CommandLineParser _commandLineParser);
  void updateFrame();
  void readFrame();
};

class InitialTimestepReader {
 public:
  InitialTimestepReader(CommandLineParser commandLineParser, AtomArray* _atomArrayPtr, SimData* _simDataPtr);
  void setContext(CommandLineParser commandLineParser, AtomArray* _atomArrayPtr, SimData* _simDataPtr);
  string initLoc;
  TimestepReader timestepReader;
  AtomArray* atomArrayPtr;
  Timestep* timestepPtr;
  InputStream inputStream;
  Timestep emptyTimestep;
  SimData* simDataPtr;
  Atom atom;

  bool fileExists();
  void openStream(CommandLineParser _commandLineParser);
  void readFromFile();
  void readFromStream();
  void writeToFile();
  void readInitialTimestep();
};

#endif
