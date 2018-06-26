#ifndef DROPLET_H
#define DROPLET_H

#include "Utils.h"
#include "Atoms.h"
#include "Fitting.h"

using namespace std;

////////////////////////
// Droplet Components //
////////////////////////

struct Monolayer
{
  double radius;
  double height;
  Options options;
  SimData *simDataPtr;
  AtomArray *atomArrayPtr;
  double zlim[2];
  TH1D *hMono;
  TanhFit tanhFit;

  Monolayer();
  ~Monolayer();
  void setContext(Options _options, SimData *_simDataPtr, AtomArray *_atomArrayPtr);
  void calculateRadius();
  void reset();
  void fillOne(Atom &atom);
  void fill(AtomArray &atoms);
  bool inMonolayer(Atom &atom);
  void convertUnits();

  //Keep track of which atoms join the monolayer, and save their radius scaled by the base radius
  int monoFlux(vector<double> r,vector<double>z,double* monoLimits,double baseRadius,TH1D* rScaledJoin,TH1D* rScaledLeave,int &nMono);

  //Find z interface between monolayer and bulk
  void findMonoLimits(TH1D *hWaterDens,double *monoLimits);

};

struct CircularBulk
{
  CircularBulk();
  ~CircularBulk();

  double firstBulkBin;
  double height;
  double radius;
  double volume;
  double contactAngle;
  CircleFit circle;
  TanhFit tanhFit;
  Options options;
  SimData *simDataPtr;
  AtomArray *atomArrayPtr;
  TGraph *gCirclePoints;
  TH2D *hDroplet;

  int numPoints = 0;
  double *boundaryPointsArray[2];
  const char* headers[2] = {"x", "y"};

  void setContext(Options options, SimData *_simDataPtr, AtomArray *_atomArrayPtr);
  void setHist(TH2D *_hDroplet);
  void fillOne(Atom &atom);
  void calculateHeight();
  void calculateRadius();
  void calculateContactAngle();
  void calculateSphericalVolume();
  void calculateCylindricalVolume();

  void saveBoundaryPoints();
  // TODO: Clean this up
  //Find boundary points by tanh fitting for each row
  void findBoundaryPoints();
  //void findBoundaryPoints(TH2D* hist,TGraph *circlePointsGraph,char* aOrR,double *monoLimits,TH1D *hMono,double& rBulkMax,double &monoEdge,double &bulkEdge,int frameStep,TLegend* tanhLegend,TLine** tanhLines,TText **tanhTexts,TPaveText *tanhTextBox);

  //Given a graph of points, a maximum x value, and the y coordinate of the interface with the substrate, fit a circle to the points
  double fitCircle(TGraph* circlePointsGraph,CircleFit &circle,double xMax,int timestep);
};

struct Droplet {
  CircularBulk bulk;
  Monolayer monolayer;
  Options options;
  SimData *simDataPtr;
  AtomArray *atomArrayPtr;
  TH2D *hDroplet;
  TH1D* hLiquidDens;
  TCanvas *cDroplet;
  double rDensCyl;

  Droplet(AtomArray &atomArray);
  ~Droplet();
  void setContext(AtomArray &atomArray);
  void fillOne(Atom &atom);
  void fill(AtomArray &atomArray);
  void convertUnits();
  void createHists();
  void createCanvas();
  void plotDensity(char* filename);
  double getMass();
  double getMass1D();
  void reset();
  void findMonolayer();
  void dropletCalculations();
};


////////////////////
// Misc. Analysis //
////////////////////

struct monolayerTracker
{
  int numMonoIDs; // # of atoms to track
  int* id; // IDs of atoms
  int* monoIDs; // ??
  AtomArray monoAtoms;
};


#endif
