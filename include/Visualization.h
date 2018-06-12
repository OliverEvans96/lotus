#ifndef VISUALIZATION_H
#define VISUALIZATION_H

#include "TLine.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TEllipse.h"
#include "TH1D.h"

#include "Utils.h"
#include "MDBase.h"
#include "Quiver.h"
#include "FieldViz.h"
#include "CircleFit.h"

using namespace std;

///////////////////
// Visualization //
///////////////////


struct Figure {
  TCanvas* canvas;
  string outFile;
  string title;
  int width;
  int height;
  SimData* simDataPtr;
  Options options;
  double xlo, xhi;
  double ylo, yhi;

  Figure(string _string, string _outFile, SimData &simData);
  ~Figure();

  void createCanvas();
  void setCanvasStyle();
  void saveImage();
  void saveROOT();
  void save();
};

struct DropletFigure : Figure {
  TH2D* hDroplet;
  TEllipse* circleEllipse;
  TGraph* gCirclePoints;
  CircleFit* circlePtr;
  TLine *tangentLine;
  TLine* bulkEdgeLine;
  TLine* monoEdgeLine;
  TLine* heightLine;
  TLine* monoHiLine;
  TLine* monoLoLine;
  TPaveText* textBox;
  TLegend* legend;
  TText* cAText;
  TText* dHText;
  TText* bEText;
  TText* mEText;

  double bulkEdge;
  double monoEdge;
  double dropletHeight;
  double contactAngle;
  double *monoLimits;

  DropletFigure(TH2D* _hDroplet, TGraph* _gCirclePoints, CircleFit &circle, string _title, string _outFile, SimData &simData);
  ~DropletFigure();

  void createLines();
  void createLegend();
  void deleteLines();
  void deleteLegend();

  void setLineStyle();
  void setHistStyle();
  void setLegendStyle();
  void setStyle();

  void setValues(double _bulkEdge, double _monoEdge, double _dropletHeight, double _contactAngle, double* _monoLimts);

  void setLegendText();
  void addLegendEntries();

  void drawHist();
  void drawLines();
  void drawText();
  void drawLegend();
  void drawAnnotations();
  void draw();
};


struct DensFigure : Figure {
  TLine* monoHiLineDens;
  TLine* monoLoLineDens;
  TH1D* hLiquidDens;
  TH1D* hSubstrateDens;
  TLegend* densLeg;
  double* monoLimits;

  DensFigure(string _string, string _outFile, SimData &simData);
  ~DensFigure();
  void setLineStyle();
  void setLegendStyle();
  void setStyle();

  void setValues(double *_monoLimits);

  void createLines();
  void createLegend();

  void deleteLines();
  void deleteLegend();

  void drawLines();
  void drawLegend();
  void draw();
};

struct TanhFigure : Figure {
  TLine* ldTanhLine;
  TLine* solTanhLine;
  TLine* x0TanhLine;
  TLine* x0guessTanhLine;
  TLine* lowGuessTanhLine;
  TLine* hiGuessTanhLine;
  TLine* lowBinTanhLine;
  TLegend* legend;
  TPaveText* textBox;
  TText* tanhTexts;
  TText* tanhLines;

  TanhFigure(string _string, string _outFile, SimData &simData);
  ~TanhFigure();
  void createLines();
  void setStyle();
  void fillLegend();
  void fillPaveText();
  void draw();
};

//Draw a TH1D horizontally, returning a TGraph which should probably be deleted
TGraph *horizontalHist(TH1D* hist);

#endif
