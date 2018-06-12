#ifndef VISUALIZATION_H
#define VISUALIZATION_H

#include "TLine.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TH1D.h"

#include "Utils.h"
#include "Quiver.h"
#include "FieldViz.h"

using namespace std;

///////////////////
// Visualization //
///////////////////


struct Figure
{
  TCanvas* canvas;
  string outFile;
  string title;
  int width;
  int height;
  SimData* simDataPtr;
  double xlo, xhi;
  double ylo, yhi;

  Figure(string _string, string _outFile);
  ~Figure();
  void save();
  void saveROOT();
};

struct DropletFigure : Figure
{
  TH2D *hDroplet;
  TEllipse *circleEllipse;
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

  DropletFigure(TH2D* _hDroplet, string _title, string _outFile, SimData* _simDataPtr);
  ~DropletFigure();

  void createLines();
  void createLegend();
  void deleteLines();
  void deleteLegend();

  void setLineStyle();
  void setHistStyle();
  void setLegendStyle();
  void setStyle();

  void drawHist();
  void drawLines();
  void drawText();
  void drawAnnotations();
  void draw();
};


class DensFigure : Figure
{
  TLine* monoHiLineDens;
  TLine* monoLoLineDens;
  TH1D* hLiquidDens;
  TH1D* hSubstrateDens;
  TLegend* densLeg;

  DensFigure();
  ~DensFigure();
  void setLineStyle();
  void setLegendStyle();
  void setStyle();
  void createLines();
  void createLegend();
  void deleteLines();
  void deleteLegend();
  void drawLines();
  void drawLegend();
  void draw();

  // TODO: Remove this
  void plotDensity();
};

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
  TPaveText* textBox;
  TText* tanhTexts;
  TText* tanhLines;

  TanhFigure();
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
