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

AtomArray::AtomArray() {
  cout << "Warning: using default constructor." << endl;
}

AtomArray::AtomArray(SimData &simData) {
  setNumAtoms(simData.numAtoms);
  allocateArrays();
}

AtomArray::~AtomArray() {
  deallocateArrays();
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

void AtomArray::deallocateArrays() {
  delete [] x;
  delete [] y;
  delete [] z;
  delete [] r;
  delete [] p;

  delete [] vx;
  delete [] vy;
  delete [] vz;
  delete [] vr;
  delete [] vp;
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

void AtomArray::getAtom(int i, Atom atom) {
  atom.x = x[i];
  atom.y = y[i];
  atom.z = z[i];
  atom.r = r[i];
  atom.p = p[i];

  atom.vx = vx[i];
  atom.vy = vy[i];
  atom.vz = vz[i];
  atom.vr = vr[i];
  atom.vp = vp[i];

  atom.cosTheta = cosTheta[i];
}

