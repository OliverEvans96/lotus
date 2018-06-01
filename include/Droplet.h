#ifndef DROPLET_H
#define DROPLET_H

#include "Utils.h"
#include "Atoms.h"
#include "CircleFit.h"

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
  double zlim[2]; //TODO: set
  // TODO: Use pointer here, initialize properly
  TH1D hMono;

  Monolayer();
  void setContext(Options _options, SimData *_simDataPtr, AtomArray *_atomArrayPtr);
  void calculateRadius();
  void fillOne(Atom &atom);
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

  double height;
  double radius;
  double volume;
  double contactAngle;
  CircleFit circle;
  Options options;
  SimData *simDataPtr;
  AtomArray *atomArrayPtr;

  void setContext(Options options, SimData *simDataPtr, AtomArray *atomArrayPtr);
  void fillOne(Atom &atom);
  void calculateHeight();
  void calculateRadius();
  void calculateContactAngle();
  void calculateSphericalVolume();
  void calculateCylindricalVolume();

  double guessRowBoundary(TH2D* hist,int j);

  //Find boundary points by tanh fitting for each row
  void findBoundaryPoints(TH2D* hist,TGraph *circlePointsGraph,char* aOrR,double *monoLimits,TH1D *hMono,double& rBulkMax,double &monoEdge,double &bulkEdge,int frameStep,TLegend* tanhLegend,TLine** tanhLines,TText **tanhTexts,TPaveText *tanhTextBox);

  //Given a graph of points, a maximum x value, and the y coordinate of the interface with the substrate, fit a circle to the points
  double fitCircle(TGraph* circlePointsGraph,CircleFit &circle,double xMax,int timestep);

  //Guess boundary of water molecule by counting for a single column
  void guessTanhFit(TH1D* hist,double* fitBounds,double &ld, double &width,double &boundary,double &xlo,double &xhi);

  double solveTanhFit(TH1D* hist, TF1* tanhFit, double* fitBounds, int startBin, int frameStep, double bulkEdge, string fitType, double pos, TLegend* tanhLegend, TLine** tanhLines, TText **tanhTexts, TPaveText *tanhTextBox);
};

struct Droplet {
  CircularBulk bulk;
  Monolayer monolayer;
  Options options;
  SimData *simDataPtr;
  AtomArray *atomArrayPtr;
  TH2D *hDroplet;
  TCanvas *cDroplet;

  Droplet(AtomArray &atomArray);
  ~Droplet();
  void setContext(AtomArray &atomArray);
  void fillOne(Atom &atom);
  void fill(AtomArray &atomArray);
  void convertUnits();
  void createHist();
  void createCanvas();
  void plotDensity(char* filename);
  double getMass();
  void reset();
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
