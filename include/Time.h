/**
   @file Time.h

   Timekeeping variables such as Timestep and Frame.

   Before any calculation occurs, multiple timesteps
   are averaged into a frame (as defined by Options::stepsPerFrame).
   The Timestep and Frame objects are declared here.

   If the number of steps in the dumpfile is not divisible
   by Options::stepsPerFrame, then the last frame contains extra
   timesteps. Relevant quantities are stored in the LastFrame object.
*/

#ifndef TIME_H
#define TIME_H

#include <iostream>

using namespace std;

//////////////////
// Time Keepers //
//////////////////

/**
   A single timestep.
*/
struct Timestep
{
  /** The timestep number from the LAMMPS dumpfile,
      corresponding to actual simulation time (generally in fs)
  */
  int time; // 1 fs
  /** The index of the timestep.
      e.g. 0, 1, 2 for the first three timesteps.
  */
  int stepNum; // 2 ps (index)

  Timestep();
};

struct Frame
{
  /// Timestep::time for the first timestep in the frame
  int time; // 1 fs (time of 1st timestep in frame)
  /** The index of the frame.
      e.g. 0, 1, 2 for the first three frames.
  */
  int frameNum; // 10 ps
  /// Timestep::stepNum for the first timestep in the frame
  int frameStep; // stepNum (index of first timestep in frame)
  /** Number of steps in the current frame.
      For all but the last frame, this is equal to
      Options::stepsPerFrame.
      For the last frame, it is equal to
      LastFrame::numSteps.
  */
  int stepsThisFrame; // same for all but last frame
  /// The number of atoms times the number of steps this frame
  int atomsThisFrame; // numAtoms * stepsThisFrame

  Frame();
};

// Deal separately with last frame if numSteps is not
// divisible by stepsPerFrame.
/**
   If Options::stepsPerFrame does not evenly divide
   SimData::numSteps, then the last frame contains extra steps.
*/
struct LastFrame
{
  /// Number of additional steps in the last frame
  int extraSteps;
  /// Frame::frameNum for the last frame
  int frameNum;
  /** Number of steps in the last frame
      @see Frame::stepsThisFrame
  */
  int numSteps;

  /// Calculate and set the number of steps in the last frame
  void setSteps(int numSteps, int stepsPerFrame);
};

#endif
