/**
   @file MDBase.h

   Basic quantities related to the MD simulation.

   Contains time-independent data such as grid and
   atom information, as well as time-dependent
   information such as the substrate and monolayer extents.
*/

#ifndef MDBASE_H
#define MDBASE_H

#include "Utils.h"
#include "Parameters.h"
#include "Time.h"

using namespace std;

/////////////////////
// Simulation Data //
/////////////////////

/**
  Simulation box bounds, as defined in the LAMMPS dumpfile.
  This data is updated each timestep.

  @see DumpfileReader
*/
struct BoxBounds {
  double xlo, xhi;
  double ylo, yhi;
  double zlo, zhi;
};

/**
  Source of truth for general MD variables.
  Variables are set here, and a pointer to this object
  is passed around and read elsewhere.
*/
struct SimData {
  // Static data
  Options options;
  int numAtoms;
  int numSteps;
  int numFrames;
  int stepsPerFrame;
  LastFrame lastFrame;
  /// Which atom types correspond to liquid
  vector<int> liquidTypes;
  /// Which atom types correspond to substrate
  vector<int> solidTypes;
  /// Mass for each atom type
  map<int, double> masses;
  /** Maps oxygen atom id to two hydrogen atom ids.
      @see DatafileReader::readWaterBonds
  */
  map<int, int*> waterBonds;
  BoxBounds simBounds;

  // Dynamic data
  /** z coordinate of the top of the substrate.
     If the substrate option is disabled,
     this can be manually specified in options.
     @todo link to actual option.
     @see Substrate::findLimits
  */
  double substrateTop;
  /** z coordinate of the top of the substrate.
     If the substrate option is disabled,
     this can be manually specified in options.
     @todo link to actual option.
  */
  double monoTop;
  Frame *framePtr;

  SimData(Options options);
  ~SimData();

  void deleteWaterBonds();
  void setOptions(Options options);
  void setNumSteps(int _numSteps);
  void setStepsPerFrame(int _stepsPerFrame);
};

struct Grid {
  double zlo, zhi, rhi, vhi;
  double dz, dv;
  int nz, nv, nr;

  double* rVals;
  double zhi_preround, vhi_preround, rhi_preround;
  bool allocated;

  ~Grid();

  void setBounds(double _zlo, double _zhi, double _rhi);
  void setSpacing(double _dz, double _dv);

  void init();
  void calculateVolumeLimits();
  void allocateBins();
  void calculateBins();
  void deallocateBins();
};

//double fitCircle(TGraph2D* circlePointsGraph,CircleFit circle,double xMax,int timestep);


//Count the number of atoms which have migrated from the bulk to the monolayer this timestep
//double bulkMonoExchange(vector<double>z,double *monoLimits,int stepsPerFrame);


//Find the lowest atom ID for a water molecule in this simulation
int lowestID(ifstream &inFile);




#endif
