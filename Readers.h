#ifndef READERS_H
#define READERS_H

using namespace std;

/////////////
// Readers //
/////////////

struct LineReader
{
  string line;
  string input;
  int atomNum;
  int lineNum;
  Atom atom;

  void readLine();
};

struct TimestepReader
{
  LineReader lineReader;
  Timestep timestep;
  int numSteps;

  void readTimestep();
};

struct FrameReader
{
  TimestepReader timestepReader;
  Frame frame;
  int stepsPerFrame;
};

#endif
