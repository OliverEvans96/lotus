#ifndef ATOMS_H
#define ATOMS_H

#include "Utils.h"

using namespace std;

//////////////////////
// Single-atom data //
//////////////////////

struct Position
{
  double x, y, z, r, p;
  void calculateNonCartesian();
};

struct Velocity
{
  double vx, vy, vz, vr, vp;
  void calculateNonCartesian();
};

struct Atom
{
  Velocity velocity;
  Position position;
};

/////////////////////
// Multi-atom data //
/////////////////////

struct PositionArray
{
  double* x, y, z, r, p;
  void calculateNonCartesian();
};

struct VelocityArray
{
  double* vx, vy, vz, vr, vp;
  void calculateNonCartesian();
};

struct AtomArray
{
  VelocityArray velocity;
  PositionArray position;
};



#endif
