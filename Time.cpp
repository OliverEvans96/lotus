#include "Time.h"

/////////////////////////
// Variables over time //
/////////////////////////

StepVariables::StepVariables()

FrameVariables::FrameVariables()


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

LastFrame::LastFrame()

LastFrame::setSteps(int numSteps, int stepsPerFrame) {
  extraSteps = numSteps % stepsPerFrame;
  divisible = (extraSteps == 0);
  penultimateFrame = numSteps - extraSteps;

  cout << "extraSteps: " << extraSteps << endl;
  cout << "divisible: " << divisible << endl;
  cout << "penultimateFrame: " << penultimateFrame << endl;

}
