#ifndef TIME_H
#define TIME_H

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
  int fs_per_step;
  StepVariables stepVariables;

  void incrementTimestep();
};

struct Frame
{
  int time; // 1 fs
  int frameNum; // 10 ps
  int frameStep; // stepNum of first timestep in frame
  bool frameNamed; // Whether 
};

// Deal separately with last frame if numSteps is not
// divisible by stepsPerFrame.
struct LastFrame
{
  int penultimateFrame;
  int extraSteps;
  bool divisible;
};


#endif
