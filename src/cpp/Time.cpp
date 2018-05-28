#include "Time.h"

/////////////////////////
// Variables over time //
/////////////////////////


//////////////////
// Time Keepers //
//////////////////

Timestep::Timestep() {
  time = 0;
  stepNum = 0;
}

Frame::Frame() {
  time = 0;
  frameNum = 0;
  frameStep = 0;
}

void LastFrame::setSteps(int totalNumSteps, int stepsPerFrame) {
  extraSteps = totalNumSteps % stepsPerFrame;
  numSteps = stepsPerFrame + extraSteps;
  frameNum = totalNumSteps / stepsPerFrame - 1;
}
