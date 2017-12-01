#ifndef MDBASE_H
#define MDBASE_H

#include <iostream>
using namespace std;

///////////////////////////
// High level parameters //
///////////////////////////

struct CommandLineParser
{
  string inLoc, outLoc;
  void parse();
}

struct Options
{
  bool skipToEnd; // Jump to last frame
  bool trackMonoAtoms; // Follow a few atoms that end up in the monolayer
  bool saveImages; // Produce images or not
  bool plotHist; // Whether to draw 2D cylindrical histogram
  bool plotDipole;
  bool plotVr;
  bool plotDensity;
  bool plotAllTogether;
  bool debugOutput; // Enable debuguging print statments
}

struct AnalysisParameters
{
  double rDensCyl; // Maximum cylinder to use for radial binning
  double rBulkMax; // ??
  double binVolume; // Constant volume for all cylindrical bins
}

///////////////////
// Visualization //
///////////////////

struct Figure
{
  TCanvas* canvas;
  string outFile;
  string title;

  void draw();
  void save();
}

struct HistFigure : Figure
{
  TLine* bulkEdgeLine;
  TLine* monoEdgeLine;
  TLine* heightLine;
  TLine* monoHiLine;
  TLine* monoLoLine;
  TTextBox* textBox;
  TLegend* legend;
  TText* cAText;
  TText* dHText;
  TText* bEText;
  TText* mEText;

  void drawText();
  void drawAnnotations();
}


class DensFigure : Figure
{
  TLine* monoHiLineDens;
  TLine* monoLoLineDens;

  void drawAnnotations();
}

class TanhFigure : Figure
{
  TLine* ldTanhLine;
  TLine* solTanhLine;
  TLine* x0TanhLine;
  TLine* x0guessTanhLine;
  TLine* lowGuessTanhLine;
  TLine* hiGuessTanhLine;
  TLine* lowBinTanhLine;
  TLegend* legend;
  TTextBox* textBox;
  TText* tanhTexts;
  TText* tanhLines;

  void setStyle();
  void fillLegend();
  void fillTextBox();
}


//////////////////////
// Single-atom data //
//////////////////////

struct Position
{
  double x, y, z, r, p;
  void calculateNonCartesian();
}

struct Velocity
{
  double vx, vy, vz, vr, vp;
  void calculateNonCartesian();
}

struct Atom
{
  Velocity velocity;
  Position position;
}

/////////////////////
// Multi-atom data //
/////////////////////

struct PositionArray
{
  double* x, y, z, r, p;
  void calculateNonCartesian();
}

struct VelocityArray
{
  double* vx, vy, vz, vr, vp;
  void calculateNonCartesian();
}

struct AtomArray
{
  VelocityArray velocity;
  PositionArray position;
}


/////////////////////////
// Variables over time //
/////////////////////////

struct StepVariables
{
  double monoEdge;
  double bulkEdge;
  double contactAngle;
  double dropletHeight;
  double avgRadius;
  double avgDipole;
  double MSD;
  double frameFlux;
  double dMe_dt;
  double nMono;
  double chi2;
}

struct FrameVariables
{
  double monoEdge;
  double bulkEdge;
  double contactAngle;
  double dropletHeight;
  double avgRadius;
  double avgDipole;
  double MSD;
  double frameFlux;
  double dMe_dt;
  double nMono;
  double chi2;
}


//////////////////
// Time Keepers //
//////////////////

struct Timestep
{
  int time; // 1 fs
  int stepNum; // 2 ps
  StepVariables stepVariables;
}

struct Frame
{
  int time; // 1 fs
  int frameNum; // 10 ps
  int frameStep; // ??
  bool frameNamed; // ??
}

struct LastFrame // ??
{
  int penultimateFrame;
  int extraSteps;
  bool divisible;
}

/////////////
// Readers //
/////////////

struct LineReader
{
  string line;
  string input;
  int atomNum;
  int lineNum;
  Atom atom;

  void readLine();
}

struct TimestepReader
{
  LineReader lineReader;
  Timestep timestep;
  int numSteps;

  void readTimestep();
}

struct FrameReader
{
  TimestepReader timestepReader;
  Frame frame;
  int stepsPerFrame;
}


////////////////////////
// Droplet Components //
////////////////////////

struct Monolayer
{
  double radius;
  double height;
  AtomArray atoms;

  void calculateRadius();
}

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
}

struct SphericalBulk : CircularBulk
{
  void calculateColume();
}

struct CylindricalBulk : CircularBulk
{
  void calculateColume();
}


////////////////////
// Misc. Analysis //
////////////////////

struct monolayerTracker
{
  int numMonoIDs; // # of atoms to track
  int* id; // IDs of atoms 
  int* monoIDs; // ??
  AtomArray monoAtoms;
}

#endif
