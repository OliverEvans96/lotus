/**
   @file Atoms.h

   @brief Structures for single and multiple atom positions.

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
   @brief Simple structure to hold information about a single atom.
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
   @brief Data structure to hold information about multiple atoms.

   Can get and set coordinates from Atom.
   This object must be instantiated after DatafileReader
   so that SimData::numAtoms is already defined
   since this value is used to determine the size of allocated arrays.

   Atom positions for all timesteps in a frame are stored,
   which may lead to memory issues if Options::stepsPerFrame
   or SimData::numAtoms are very large.
*/
class AtomArray
{
  /// @brief Whether the position arrays have been allocated
  bool allocated;
  void deallocateArrays();
  /// @brief Return index of atom position in arrays.
  int getIndex(int atomNum, int stepInFrame);

 public:
  void allocateArrays();
  int numAtoms;
  SimData *simDataPtr;
  int *type;
  double *x, *y, *z, *r;

  AtomArray(SimData &simData);
  ~AtomArray();

  /// @brief Set simDataPtr and number of atoms.
  void setSimData(SimData &simData);
  void setNumAtoms(int _numAtoms);
  /// @brief Copy data for atom @p atomNum at step @p stepInFrame to @p atom.
  void setAtom(int atomNum, int stepInFrame, Atom &atom);
  /// @brief Copy data for atom @p atomNum at step @p stepInFrame from @p atom.
  void getAtom(int atomNum, int stepInFrame, Atom &atom);

  /// @brief Print min, max, mean, std. of each x, y, z components of positions.
  void printStats();
};

#endif
