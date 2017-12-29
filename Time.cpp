#include "Time.h"

/////////////////////////
// Variables over time //
/////////////////////////

StepVariables::StepVariables()

FrameVariables::FrameVariables()


//////////////////
// Time Keepers //
//////////////////

Timestep::Timestep(time, fs_per_step)
{
  this->time = time;
  this->fs_per_step = fs_per_step
    }

Frame::Frame()

LastFrame::LastFrame ()

LastFrame::setSteps(int numSteps, int stepsPerFrame) {
  extraSteps = numSteps % stepsPerFrame;
  divisible = (extraSteps == 0);
  penultimateFrame = numSteps - extraSteps;
}
