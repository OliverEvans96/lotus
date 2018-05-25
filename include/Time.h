#ifndef TIME_H
#define TIME_H

#include <iostream>

using namespace std;

/////////////////////////
// Variables over time //
/////////////////////////

struct StepVariables
{
  double monoEdge;
  double bulkEdge;
  double contactAngle;
  double dropletHeight;
  double avgRadius;
  double avgDipole;
  double MSD;
  double frameFlux;
  double dMe_dt;
  double nMono;
  double chi2;
};

struct FrameVariables
{
  double monoEdge;
  double bulkEdge;
  double contactAngle;
  double dropletHeight;
  double avgRadius;
  double avgDipole;
  double MSD;
  double frameFlux;
  double dMe_dt;
  double nMono;
  double chi2;
};


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

  Frame();
};

// Deal separately with last frame if numSteps is not
// divisible by stepsPerFrame.
struct LastFrame
{
  int penultimateFrame;
  int extraSteps;
  bool divisible;

  void setSteps(int numSteps, int stepsPerFrame);
};

#endif
