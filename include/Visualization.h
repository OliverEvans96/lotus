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

  void draw();
  void save();
};

struct HistFigure : Figure
{
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

  void drawText();
  void drawAnnotations();
};


class DensFigure : Figure
{
  TLine* monoHiLineDens;
  TLine* monoLoLineDens;

  void drawAnnotations();
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

  void setStyle();
  void fillLegend();
  void fillPaveText();
};

//Draw a TH1D horizontally, returning a TGraph which should probably be deleted
TGraph *horizontalHist(TH1D* hist);

#endif
