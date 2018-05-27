#ifndef ATOMS_H
#define ATOMS_H

#include <cstring>
#include <cstdio>
#include "Utils.h"
#include "MDBase.h"

using namespace std;

//////////////////////
// Single-atom data //
//////////////////////

struct Atom
{
 public:
  double x, y, z, r, p;
  double vx, vy, vz, vr, vp;
  double cosTheta;
  void calculateNonCartesian();
};

/////////////////////
// Multi-atom data //
/////////////////////

class AtomArray
{
  SimData *simDataPtr;
  void allocateArrays();
  void deallocateArrays();

 public:
  int numAtoms;
  double *x, *y, *z, *r, *p;
  double *vx, *vy, *vz, *vr, *vp;
  double *cosTheta;

  AtomArray();
  AtomArray(SimData &simData);
  ~AtomArray();
  void setSimData(SimData &simData);
  void setAtom(int i, Atom atom);
  void getAtom(int i, Atom atom);
  void setNumAtoms(int _numAtoms);
};

#endif
