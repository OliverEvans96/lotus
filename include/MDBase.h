#ifndef MDBASE_H
#define MDBASE_H

#include "Utils.h"
#include "Parameters.h"
#include "Time.h"

using namespace std;

// Conversion factor - number density to water mass density
const double convFact = 18/.60221409;

/////////////////////
// Simulation Data //
/////////////////////

struct BoxBounds {
  double xlo, xhi;
  double ylo, yhi;
  double zlo, zhi;
};

// Source of truth for general MD variables.
// Variables are set here, and this object (or a pointer)
// is passed around and read elsewhere.
struct SimData {
  // Static data
  Options options;
  int numAtoms;
  int numSteps;
  int numFrames;
  int stepsPerFrame;
  LastFrame lastFrame;
  vector<int> liquidTypes;
  vector<int> solidTypes;
  map<int, double> masses;
  map<int, int*> waterBonds;
  BoxBounds simBounds;

  // Dynamic data
  double substrateTop;
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
