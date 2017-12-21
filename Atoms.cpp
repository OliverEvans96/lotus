#include "Atoms.h"

//////////////////////
// Single-atom data //
//////////////////////

void Atom::calculateNonCartesian() {
  r = sqrt(x*x + y*y);
  p = sqrt(x*x + y*y + z*z);

  // SHOULD DOUBLE CHECK THESE - NOT 100% SURE
  vr = (vx*x + vy*y) / r;
  vp = (vx*x + vy*y + vz*z) / p;
}

/////////////////////
// Multi-atom data //
/////////////////////

AtomArray::~AtomArray() {
  delete x, y, z, r, p;
  delete vx, vy, vz, vr, vp;
}

void AtomArray::allocateArrays() {
  x = new double[numAtoms];
  y = new double[numAtoms];
  z = new double[numAtoms];
  r = new double[numAtoms];
  p = new double[numAtoms];

  vx = new double[numAtoms];
  vy = new double[numAtoms];
  vz = new double[numAtoms];
  vr = new double[numAtoms];
  vp = new double[numAtoms];
}

void AtomArray::setNumAtoms(int _numAtoms) {
  numAtoms = _numAtoms;
  allocateArrays();
}

void AtomArray::setAtom(int i, Atom atom) {
  x[i] = atom.x;
  y[i] = atom.y;
  z[i] = atom.z;
  r[i] = atom.r;
  p[i] = atom.p;

  vx[i] = atom.vx;
  vy[i] = atom.vy;
  vz[i] = atom.vz;
  vr[i] = atom.vr;
  vp[i] = atom.vp;

  cosTheta[i] = atom.cosTheta;
}


