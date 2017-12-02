#ifndef DROPLET_H
#define DROPLET_H

#include <iostream>

using namespace std;

////////////////////////
// Droplet Components //
////////////////////////

struct Monolayer
{
  double radius;
  double height;
  AtomArray atoms;

  void calculateRadius();
};

struct CircularBulk
{
  double height;
  double radius;
  double volume;
  double contactAngle;
  CircleFitClass circle;

  void calculateHeight();
  void calculateRadius();
  void calculateContactAngle();
};

struct SphericalBulk : CircularBulk
{
  void calculateColume();
};

struct CylindricalBulk : CircularBulk
{
  void calculateColume();
};

////////////////////
// Misc. Analysis //
////////////////////

struct monolayerTracker
{
  int numMonoIDs; // # of atoms to track
  int* id; // IDs of atoms 
  int* monoIDs; // ??
  AtomArray monoAtoms;
};


#endif
