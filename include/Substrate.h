#ifndef SUBSTRATE_H
#define SUBSTRATE_H

#include "Utils.h"
#include "Atoms.h"

using namespace std;

struct Substrate {
  // TODO: Use pointer here, initialize properly
  TH1D *hSubs;
  TCanvas *cSubs;
  Options options;
  SimData *simDataPtr;
  AtomArray *atomArrayPtr;
  double zlim[2]; //TODO: Set

  Substrate(AtomArray &atomArray, double dz=0.5);
  ~Substrate();
  void setContext(AtomArray &atomArray);
  void fillOne(Atom &atom);
  void fill(AtomArray &atoms);
  void reset();
  void convertUnits();
  void createHist(double dz);
  void createCanvas();
  void plotDensity(char* filename);
  void findLimits();
  double getMass();
};

#endif
