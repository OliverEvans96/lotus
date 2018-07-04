/**
   @file Time.h

   @brief Timekeeping variables such as Timestep and Frame.

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
   @brief A single timestep.
*/
struct Timestep
{
  /** @brief The timestep number from the LAMMPS dumpfile,
      corresponding to actual simulation time (generally in fs).
  */
  int time;
  /** @brief The index of the timestep.

      e.g. 0, 1, 2 for the first three timesteps.
  */
  int stepNum;

  Timestep();
};

struct Frame
{
  /// @brief Timestep::time for the first timestep in the frame
  int time;
  /** @brief The index of the frame.

      e.g. 0, 1, 2 for the first three frames.
  */
  int frameNum;
  /// @brief Timestep::stepNum for the first timestep in the frame
  int frameStep;
  /** @brief Number of steps in the current frame.

      For all but the last frame, this is equal to
      Options::stepsPerFrame.
      For the last frame, it is equal to
      LastFrame::numSteps.
  */
  int stepsThisFrame;
  /** @brief The number of atoms times the number of steps this frame.

      Equal to `numAtoms * stepsThisFrame`
  */
  int atomsThisFrame;

  Frame();
};

/**
   @brief If Options::stepsPerFrame does not evenly divide
   SimData::numSteps, then the last frame contains extra steps.
*/
struct LastFrame
{
  /// @brief Number of additional steps in the last frame
  int extraSteps;
  /// @brief Frame::frameNum for the last frame
  int frameNum;
  /** @brief Number of steps in the last frame
      @see Frame::stepsThisFrame
  */
  int numSteps;

  /// @brief Calculate and set the number of steps in the last frame
  void setSteps(int numSteps, int stepsPerFrame);
};

#endif
