#ifndef READERS_H
#define READERS_H

#include "Utils.h"
#include "Atoms.h"
#include "Time.h"

using namespace std;

struct InputStream {
  InputStream(char* _filename);
  ~InputStream();
  void verifyStream();

 public:
  char* filename;
  ifstream stream;
};


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
  void setContext(InputStream* _inputStreamPtr, int* _atomNumPtr, int* _lineNumPtr);
  void readLine();
};

class TimestepReader {
  LineReader lineReader;
  HeaderReader headerReader;
  Timestep timestep;
  AtomArray* atomArrayPtr;
  InputStream* inputStreamPtr;
  int atomNum;
  int lineNum;

 public:
  void setContext(InputStream* _inputStreamPtr, AtomArray* _atomArrayPtr);
  void resetAtomCounter();
  void readTimestep();
};

class FrameReader {
  FrameReader(InputStream* _inputStreamPtr, AtomArray* _atomArrayPtr);
  void setContext(InputStream* _inputStreamPtr, AtomArray* _atomArrayPtr);

  TimestepReader timestepReader;
  Frame frame;
  InputStream* inputStreamPtr;
  AtomArray* atomArrayPtr;

 public:
  int stepsPerFrame;
  void setStepsPerFrame(int _stepsPerFrame);
  void readFrame();
};

#endif
