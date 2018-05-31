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

void Atom::print() {
  printf("Atom: type %d @ (%.2f, %.2f, %.2f)\n", type, x, y, z);
}

/////////////////////
// Multi-atom data //
/////////////////////

AtomArray::AtomArray(SimData &simData) {
  allocated = false;
  setSimData(simData);
}

AtomArray::~AtomArray() {
  deallocateArrays();
}

void AtomArray::setSimData(SimData &simData) {
  simDataPtr = &simData;
  setNumAtoms(simDataPtr->numAtoms);
  cout << "numAtoms = " << numAtoms << endl;
}
void AtomArray::allocateArrays() {
  // Last frame may be longer than the rest,
  // so allocate arrays accordingly.
  int arraySize = numAtoms * simDataPtr->lastFrame.numSteps;
  if(!allocated ) {
    if(simDataPtr->options.verbose) {
      cout << "lastFrame.numSteps = " << simDataPtr->lastFrame.numSteps << endl;
      cout << "Allocating arrays of size " << arraySize;
      printf(" (%8.4f) GB\n", arraySize/pow(2,27)*3);
    }


    type = new int[arraySize];
    x = new double[arraySize];
    y = new double[arraySize];
    z = new double[arraySize];
    r = new double[arraySize];
    p = new double[arraySize];

    // vx = new double[numAtoms];
    // vy = new double[numAtoms];
    // vz = new double[numAtoms];
    // vr = new double[numAtoms];
    // vp = new double[numAtoms];

    allocated = true;
  }
}

void AtomArray::deallocateArrays() {
  if(allocated) {

    delete [] type;

    delete [] x;
    delete [] y;
    delete [] z;
    delete [] r;
    delete [] p;

    // delete [] vx;
    // delete [] vy;
    // delete [] vz;
    // delete [] vr;
    // delete [] vp;
  }
}

void AtomArray::setNumAtoms(int _numAtoms) {
  numAtoms = _numAtoms;
}

void AtomArray::setAtom(int i, Atom &atom) {
  type[i] = atom.type;

  x[i] = atom.x;
  y[i] = atom.y;
  z[i] = atom.z;
  r[i] = atom.r;
  p[i] = atom.p;

  // vx[i] = atom.vx;
  // vy[i] = atom.vy;
  // vz[i] = atom.vz;
  // vr[i] = atom.vr;
  // vp[i] = atom.vp;
}

void AtomArray::getAtom(int i, Atom &atom) {
  atom.type = type[i];

  atom.x = x[i];
  atom.y = y[i];
  atom.z = z[i];
  atom.r = r[i];
  atom.p = p[i];

  // atom.vx = vx[i];
  // atom.vy = vy[i];
  // atom.vz = vz[i];
  // atom.vr = vr[i];
  // atom.vp = vp[i];
}
