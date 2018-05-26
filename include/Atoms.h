#ifndef ATOMS_H
#define ATOMS_H

#include <cstring>
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
  void allocateArrays();
  void deallocateArrays();

 public:
  AtomArray();
  AtomArray(SimData &simData);
  ~AtomArray();

  int numAtoms;
  double *x, *y, *z, *r, *p;
  double *vx, *vy, *vz, *vr, *vp;
  double *cosTheta;
  void setAtom(int i, Atom atom);
  void getAtom(int i, Atom atom);
  void setNumAtoms(int _numAtoms);
};

#endif
