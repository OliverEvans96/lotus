#ifndef VISUALIZATION_H
#define VISUALIZATION_H

#include "TLine.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TEllipse.h"
#include "TH1D.h"

#include "Utils.h"
#include "Substrate.h"
#include "Droplet.h"
#include "MDBase.h"
#include "Quiver.h"
#include "FieldViz.h"
#include "Fitting.h"

using namespace std;

///////////////////
// Visualization //
///////////////////


class Figure {
 protected:
  TCanvas* canvas;
  TLegend* legend;
  string title;
  int width;
  int height;
  SimData* simDataPtr;
  Options options;
  double xlo, xhi;
  double ylo, yhi;

 public:
  Figure(string _title, SimData &simData);
  ~Figure();

  void createCanvas();
  void setCanvasStyle();
  void saveImage();
  void saveROOT();
  void save(char* filename);
};

class DropletFigure : public Figure {
  Droplet* dropletPtr;

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
  TText* cAText;
  TText* dHText;
  TText* bEText;
  TText* mEText;

  double bulkEdge;
  double monoEdge;
  double dropletHeight;
  double contactAngle;
  double *monoLimits;

 public:
  DropletFigure(string _title, Droplet &droplet);
  ~DropletFigure();

  void createLines();
  void createLegend();
  void deleteLines();
  void deleteLegend();

  void setLineStyle();
  void setHistStyle();
  void setLegendStyle();
  void setStyle();

  void setValues();

  void setLegendText();
  void addLegendEntries();

  void drawHist();
  void drawLines();
  void drawText();
  void drawLegend();
  void drawAnnotations();
  void draw();
};


class DensFigure : public Figure {
  TLine* monoHiLineDens;
  TLine* monoLoLineDens;
  TH1D* hLiquidDens;
  TH1D* hSubstrateDens;
  double monoLimits[2];

  Droplet* dropletPtr;
  Substrate* substratePtr;

 public:
  // TODO: Set title later
  DensFigure(string _title, Droplet &droplet, Substrate &substrate);
  ~DensFigure();
  void setLineStyle();
  void setLegendStyle();
  void setStyle();

  void setValues();

  void createLines();
  void createLegend();

  void deleteLines();
  void deleteLegend();

  void drawHists();
  void drawLines();
  void drawLegend();
  void draw();
};

class TanhFigure : public Figure {
  TGraph* gPoints;

  TLine* ldLine;
  TLine* halfLdLine;
  TLine* x0Line;

  TPaveText* tanhTextBox;
  TText* posText;

  double ld;
  double w_guess;
  double w;
  double x0;
  double x0_guess;
  double lowGuess;
  double hiGuess;
  double startPoint;
  double val;

  static const int numLines = 3;
  static const int numTexts = 5;

  TText* tanhTexts[numTexts];
  TLine* tanhLines[numLines] = {
    ldLine,
    halfLdLine,
    x0Line,
  };

  TanhFit* tanhFitPtr;

 public:
  TanhFigure(string _string, TanhFit &tanhFit);
  ~TanhFigure();

  void createLines();
  void createGraph();
  void createLegend();
  void deleteLines();
  void deleteGraph();
  void deleteLegend();

  void addLegendEntries();
  void setValues();
  void setText();
  void setLinePositions();

  void setHistStyle();
  void setLineStyle();
  void setGraphStyle();
  void setLegendStyle();
  void setStyle();

  void drawHist();
  void drawFunction();
  void drawLines();
  void drawGraph();
  void drawLegend();
  void draw();
};

#endif
