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
  int type;
  void calculateNonCartesian();
  void print();
};

/////////////////////
// Multi-atom data //
/////////////////////

class AtomArray
{
  bool allocated;
  void deallocateArrays();

 public:
  void allocateArrays();
  int numAtoms;
  SimData *simDataPtr;
  int *type;
  double *x, *y, *z, *r, *p;
  // double *vx, *vy, *vz, *vr, *vp;

  AtomArray(SimData &simData);
  ~AtomArray();
  void setSimData(SimData &simData);
  void setNumAtoms(int _numAtoms);
  int getIndex(int atomNum, int stepInFrame);
  void setAtom(int atomNum, int stepInFrame, Atom &atom);
  void getAtom(int atomNum, int stepInFrame, Atom &atom);
  void printStats();
};

#endif
