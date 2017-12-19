#ifndef MDBASE_H
#define MDBASE_H

#include "Utils.h"

using namespace std;

const double PI = 3.141592653589793;
//Conversion factor - number density to water mass density
const double convFact=18/.60221409;


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
