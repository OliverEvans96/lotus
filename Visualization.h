#ifndef VISUALIZATION_H
#define VISUALIZATION_H

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
  TTextBox* textBox;
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
  TTextBox* textBox;
  TText* tanhTexts;
  TText* tanhLines;

  void setStyle();
  void fillLegend();
  void fillTextBox();
};



#endif
