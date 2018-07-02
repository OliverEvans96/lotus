/**
   @file Atoms.h

   Structures for single and multiple atom positions.

   These objects facilitate storage and manipulation of atom positions
   (and possibly velocities, but not currently), and calculation of
   non-cartesian components of stored quantities.
*/

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

/**
   Simple structure to hold information about a single atom.
*/
struct Atom
{
 public:
  double x, y, z, r;
  int type;
  void calculateRadius();
  void print();
};

/////////////////////
// Multi-atom data //
/////////////////////

/**
   Data structure to hold information about multiple atoms.

   Can get and set coordinates from Atom.
   This object must be instantiated after DatafileReader
   so that `simDataPtr->numAtoms` is already defined
   since this value is used to determine the size of allocated arrays.

   Atom positions for all timesteps in a frame are stored,
   which may lead to memory issues if options.stepsPerFrame
   or `simDataPtr->numAtoms` are very large.
*/
class AtomArray
{
  /// Whether the position arrays have been allocated
  bool allocated;
  void deallocateArrays();
  /// Return index of atom position in arrays.
  int getIndex(int atomNum, int stepInFrame);

 public:
  void allocateArrays();
  int numAtoms;
  SimData *simDataPtr;
  int *type;
  double *x, *y, *z, *r;

  AtomArray(SimData &simData);
  ~AtomArray();

  /// Set simDataPtr and number of atoms.
  void setSimData(SimData &simData);
  void setNumAtoms(int _numAtoms);
  /// Copy data for atom `atomNum` at step `stepInFrame` to @p atom.
  void setAtom(int atomNum, int stepInFrame, Atom &atom);
  /// Copy data for atom `atomNum` at step `stepInFrame` from @p atom.
  void getAtom(int atomNum, int stepInFrame, Atom &atom);

  /// Print min, max, mean, std. of each x, y, z components of positions.
  void printStats();
};

#endif
