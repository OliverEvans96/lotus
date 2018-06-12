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
  TLegend* legend;
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
  TGraph* gPoints;

  TLine* ldLine;
  TLine* solLine;
  TLine* x0Line;
  TLine* x0guessLine;
  TLine* lowGuessLine;
  TLine* hiGuessLine;
  TLine* lowBinLine;

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

  TText* tanhTexts[5];
  TLine* tanhLines[7] = {
    ldLine,
    solLine,
    x0Line,
    x0guessLine,
    lowGuessLine,
    hiGuessLine,
    lowBinLine
  };

  TanhFigure(string _string, string _outFile, SimData &simData);
  ~TanhFigure();

  void createLines();
  void createGraph();
  void createLegend();
  void deleteLines();
  void deleteGraph();
  void deleteLegend();

  void addLegendEntries();
  void setValues(double _ld, double _width, double _w, double _boundary, double _x0);
  void setText();
  void setLinePosition();

  void setLineStyle();
  void setPointStyle();
  void setLegendStyle();
  void setStyle();

  void drawLines();
  void drawPoints();
  void drawLegend();
  void draw();
};


// TODO: What is this?
//Draw a TH1D horizontally, returning a TGraph which should probably be deleted
TGraph *horizontalHist(TH1D* hist);

#endif
