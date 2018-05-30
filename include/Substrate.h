#ifndef SUBSTRATE_H
#define SUBSTRATE_H

#include "Utils.h"
#include "Atoms.h"

using namespace std;

struct Substrate {
  TH1D hSubs;
  Options options;
  SimData *simDataPtr;
  AtomArray *atomArrayPtr;
  double zLim[2]; //TODO: Set

  Substrate(AtomArray &atomArray);
  void setContext(AtomArray &atomArray);
};

#endif
