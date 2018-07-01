/**
   @file    Atoms.h
   @author  Oliver Evans
   @date    7/1/2018
   @version 1.0

   @brief Structures for atom data.

   @section DESCRIPTION
   Structures which facilitate storage and manipulation of atom positions
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
   @brief Simple structure to hold information about a single atom
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
   @brief Simple structure to hold information about multiple atoms.
*/
class AtomArray
{
  bool allocated;
  void deallocateArrays();

 public:
  void allocateArrays();
  int numAtoms;
  SimData *simDataPtr;
  int *type;
  double *x, *y, *z, *r;

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
