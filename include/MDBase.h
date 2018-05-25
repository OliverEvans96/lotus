#ifndef MDBASE_H
#define MDBASE_H

#include "Utils.h"
#include "Parameters.h"
#include "Time.h"

using namespace std;

//Conversion factor - number density to water mass density
const double convFact=18/.60221409;

// Source of truth for general MD variables.
// Variables are set here, and this object (or a pointer)
// is passed around and read elsewhere.
struct SimData {
  int numAtoms;
  int numSteps;
  int numFrames;
  int stepsPerFrame;
  InputStream inputStream;
  LastFrame lastFrame;

  SimData(Options options);
  int countAtoms(ifstream &inFile);
  int countSteps(ifstream &inFile);

  void setStepsPerFrame(int stepsPerFrame);
};


struct Grid {
  double zlo, zhi, rhi, vhi;
  double dz, dv;
  int nz, nv, nr;

  double* rVals;

  double zhi_preround, vhi_preround, rhi_preround;

  ~Grid();

  void setBounds(double _zlo, double _zhi, double _rhi);
  void setSpacing(double _dz, double _dv);

  void calculateVolumeLimits();
  void allocateBins();
  void calculateBins();
  void init();
  void deallocateBins();
};

//Count the number of water atoms in the first timestep
int countAtoms(ifstream &inFile);


//double fitCircle(TGraph2D* circlePointsGraph,CircleFit circle,double xMax,int timestep);


//Calculate MSD - mean squared displacement
vector<double> calculateMSD(vector<double> xi,vector<double> yi,vector<double> zi,vector<double> x, vector<double> y,vector<double> z);

//Count the number of atoms which have migrated from the bulk to the monolayer this timestep
//double bulkMonoExchange(vector<double>z,double *monoLimits,int stepsPerFrame);


//Find the lowest atom ID for a water molecule in this simulation
int lowestID(ifstream &inFile);




#endif
