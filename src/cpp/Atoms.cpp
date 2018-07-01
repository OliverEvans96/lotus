#include "Atoms.h"

//////////////////////
// Single-atom data //
//////////////////////

/**
   @brief Calculate cylindrical radius from x and y coordinates.
*/
void Atom::calculateRadius() {
  r = sqrt(x*x + y*y);
}

/**
   @brief Print atom coordinates
*/
void Atom::print() {
  printf("Atom: type %d @ (%.2f, %.2f, %.2f)\n", type, x, y, z);
}

/////////////////////
// Multi-atom data //
/////////////////////

// Must be created after DatafileReader
// so that numAtoms is known
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
      printf(" (%8.4f GB)\n", arraySize/pow(2.0,27)*3);
    }


    type = new int[arraySize];
    x = new double[arraySize];
    y = new double[arraySize];
    z = new double[arraySize];
    r = new double[arraySize];

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
  }
}

void AtomArray::setNumAtoms(int _numAtoms) {
  numAtoms = _numAtoms;
}

int AtomArray::getIndex(int atomNum, int stepInFrame) {
  return stepInFrame * numAtoms + atomNum;
}

void AtomArray::setAtom(int atomNum, int stepInFrame, Atom &atom) {
  int i = getIndex(atomNum, stepInFrame);
  type[i] = atom.type;

  x[i] = atom.x;
  y[i] = atom.y;
  z[i] = atom.z;
  r[i] = atom.r;
}

void AtomArray::getAtom(int atomNum, int stepInFrame, Atom &atom) {
  int i = getIndex(atomNum, stepInFrame);
  atom.type = type[i];

  atom.x = x[i];
  atom.y = y[i];
  atom.z = z[i];
  atom.r = r[i];
}

void AtomArray::printStats() {
  int N = simDataPtr->framePtr->atomsThisFrame;
  cout << "===ATOM STATS===" << endl;
  printf("min: %.2f, %.2f, %.2f\n", min(x,N), min(y,N), min(z,N));
  printf("max: %.2f, %.2f, %.2f\n", max(x,N), max(y,N), max(z,N));
  printf("mean: %.2f, %.2f, %.2f\n", mean(x,N), mean(y,N), mean(z,N));
  printf("stddev: %.2f, %.2f, %.2f\n", stddev(x,N), stddev(y,N), stddev(z,N));
  cout << "================" << endl;
}
