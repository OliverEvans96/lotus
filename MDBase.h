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

  SimData(CommandLineParser commandLineParser);
  int countAtoms(ifstream &inFile);
  int countSteps(ifstream &inFile);

  void setStepsPerFrame(int stepsPerFrame);
};

//Count the number of water atoms in the first timestep
int countAtoms(ifstream &inFile);


//void findBoundaryPoints(TH2D* hist,TGraph2D *circlePointsGraph,char* aOrR,double *monoLimits,TH1D *hMono,double& rBulkMax,double &monoEdge,int frameStep);

//double fitCircle(TGraph2D* circlePointsGraph,CircleFit circle,double xMax,int timestep);


//Calculate MSD - mean squared displacement
vector<double> calculateMSD(vector<double> xi,vector<double> yi,vector<double> zi,vector<double> x, vector<double> y,vector<double> z);

//Count the number of atoms which have migrated from the bulk to the monolayer this timestep
//double bulkMonoExchange(vector<double>z,double *monoLimits,int stepsPerFrame);


//Find the lowest atom ID for a water molecule in this simulation
int lowestID(ifstream &inFile);




#endif
