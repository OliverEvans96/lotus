#ifndef TIME_H
#define TIME_H

#include <iostream>

using namespace std;

//////////////////
// Time Keepers //
//////////////////

struct Timestep
{
  int time; // 1 fs
  int stepNum; // 2 ps (index)
  StepVariables stepVariables;

  Timestep();

  void incrementTimestep();
};

struct Frame
{
  int time; // 1 fs (time of 1st timestep in frame)
  int frameNum; // 10 ps
  int frameStep; // stepNum (index of first timestep in frame)
  int stepsThisFrame; // same for all but last frame
  int atomsThisFrame; // numAtoms * stepsThisFrame

  Frame();
};

// Deal separately with last frame if numSteps is not
// divisible by stepsPerFrame.
struct LastFrame
{
  int extraSteps;
  int frameNum;
  int numSteps;

  void setSteps(int numSteps, int stepsPerFrame);
};

#endif
