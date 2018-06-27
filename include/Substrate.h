#ifndef SUBSTRATE_H
#define SUBSTRATE_H

#include "Utils.h"
#include "Atoms.h"

using namespace std;

struct Substrate {
  TH1D *hSubstrateDens;
  Options options;
  SimData *simDataPtr;
  AtomArray *atomArrayPtr;
  double zlim[2]; //TODO: Set
  double dz;
  double rDensCyl;

  Substrate(AtomArray &atomArray, double dz=0.5);
  ~Substrate();
  void setContext(AtomArray &atomArray);
  void fillOne(Atom &atom);
  void fill(AtomArray &atoms);
  void reset();
  void convertUnits();
  void createHist();
  void findLimits();
  double getMass();
};

#endif
