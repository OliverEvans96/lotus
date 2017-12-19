#ifndef DROPLET_H
#define DROPLET_H

#include "Utils.h"
#include "CircleFitClass.h"

using namespace std;

////////////////////////
// Droplet Components //
////////////////////////

struct Monolayer
{
  double radius;
  double height;
  AtomArray atoms;

  void calculateRadius();

  //Keep track of which atoms join the monolayer, and save their radius scaled by the base radius
  int monoFlux(vector<double> r,vector<double>z,double* monoLimits,double baseRadius,TH1D* rScaledJoin,TH1D* rScaledLeave,int &nMono);

  //Find z interface between monolayer and bulk
  void findMonoLimits(TH1D *hWaterDens,double *monoLimits);

};

struct CircularBulk
{
  double height;
  double radius;
  double volume;
  double contactAngle;
  CircleFitClass circle;

  void calculateHeight();
  void calculateRadius();
  void calculateContactAngle();

  double guessRowBOundary(TH2D* hist,int j);

  //Find boundary points by tanh fitting for each row
  void findBoundaryPoints(TH2D* hist,TGraph *circlePointsGraph,char* aOrR,double *monoLimits,TH1D *hMono,double& rBulkMax,double &monoEdge,double &bulkEdge,int frameStep,TLegend* tanhLegend,TLine** tanhLines,TText **tanhTexts,TPaveText *tanhTextBox);

  //Given a graph of points, a maximum x value, and the y coordinate of the interface with the substrate, fit a circle to the points
  double fitCircle(TGraph* circlePointsGraph,CircleFit &circle,double xMax,int timestep);

  //Guess boundary of water molecule by counting for a single column
  void guessTanhFit(TH1D* hist,double* fitBounds,double &width,double &boundary,double &xlo,double &xhi);

  double solveTanhFit(TH1D* hist,TF1* tanhFit,double* fitBounds,int startBin,int frameStep,double bulkEdge,char *fitType,double pos,TLegend* tanhLegend,TLine** tanhLines,TText **tanhTexts,TPaveText *tanhTextBox);

};

struct SphericalBulk : CircularBulk
{
  void calculateColume();
};

struct CylindricalBulk : CircularBulk
{
  void calculateColume();
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
