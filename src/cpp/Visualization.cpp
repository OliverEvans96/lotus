
#include "Visualization.h"

///////////////////
// Visualization //
///////////////////

Figure::Figure(string _title, string _outFile, SimData &simData) {
  title = _title;
  outFile = _outFile;
  simDataPtr = &simData;
  options = simDataPtr->options;

  createCanvas();
  setCanvasStyle();
}

Figure::~Figure() {
  delete canvas;
}

void Figure::createCanvas() {
  width = options.plot_width;
  height = (int) round(width / options.plot_aspect);
  canvas = new TCanvas(title.data(), title.data(), width, height);
}

void Figure::setCanvasStyle() {
  // Use OpenGL for antialiasing
  gStyle->SetCanvasPreferGL(true);
}

void Figure::saveImage() {
  canvas->SaveAs(outFile.data());
}

void Figure::saveROOT() {
  // TODO: save w/ .C extension
}

// TODO: Take filename argument
// TODO: Save relative to image directory
void Figure::save(char* filename) {
  canvas->SaveAs(filename);
  // if(options.saveImages) {
  //   saveImage();
  // }

  // if(options.saveROOT) {
  //   saveROOT();
  // }
}

DropletFigure::DropletFigure(string _title, string _outFile, Droplet &droplet) : Figure(_title, _outFile, *droplet.simDataPtr){
  dropletPtr = &droplet;
  hDroplet = dropletPtr->hDroplet;
  gCirclePoints = dropletPtr->bulk.gCirclePoints;
  circlePtr = &dropletPtr->bulk.circle;

  xlo = 0;
  xhi = options.plot_rmax;
  ylo = 0;
  yhi = options.plot_zmax;

  createLines();
  createLegend();
}

DropletFigure::~DropletFigure() {
  deleteLines();
  deleteLegend();
}

void DropletFigure::createLines() {
  //Bulk & Mono radius vertical lines
  bulkEdgeLine = new TLine(xlo,ylo,xlo,yhi); //Bulk edge
  monoEdgeLine = new TLine(xlo,ylo,xlo,yhi); //Mono edge
  heightLine = new TLine(xlo,ylo,xlo,yhi); //Droplet height
  monoHiLine = new TLine(xlo,ylo,xlo,yhi); //top of monolayer
  monoLoLine = new TLine(xlo,ylo,xlo,yhi); //bottom of monolayer
}

void DropletFigure::createLegend() {
  legend = new TLegend(.65,.65,.85,.85);
  textBox = new TPaveText();
  cAText = textBox->AddText("Contact angle");
  dHText = textBox->AddText("Droplet height");
  bEText = textBox->AddText("Bulk edge");
  mEText = textBox->AddText("Mono edge");
}

void DropletFigure::deleteLines() {
  delete bulkEdgeLine;
  delete monoEdgeLine;
  delete heightLine;
  delete monoHiLine;
  delete monoLoLine;
}

void DropletFigure::deleteLegend() {
  delete legend;
  delete textBox;
  delete cAText;
  delete dHText;
  delete bEText;
  delete mEText;
}

void DropletFigure::setLineStyle() {
  //Line properties
  bulkEdgeLine->SetLineWidth(3);
  bulkEdgeLine->SetLineColor(kGreen);
  monoEdgeLine->SetLineWidth(3);
  monoEdgeLine->SetLineColor(kRed);

  heightLine->SetLineWidth(3);
  heightLine->SetLineColor(kOrange+3);

  monoHiLine->SetLineWidth(3);
  monoHiLine->SetLineColor(kBlue);
  monoHiLine->SetLineStyle(9);
  monoLoLine->SetLineWidth(3);
  monoLoLine->SetLineColor(kBlue);
  monoLoLine->SetLineStyle(2);
}

void DropletFigure::setHistStyle() {
  // TODO: Set limits from options
  double colzMin=0.0;
  double colzMax=1.5;
  hDroplet->SetMinimum(colzMin);
  hDroplet->SetMaximum(colzMax);
}

void DropletFigure::setLegendStyle() {
  //Text to show data
  textBox->SetX1NDC(.65);
  textBox->SetY1NDC(.5);
  textBox->SetX2NDC(.85);
  textBox->SetY2NDC(.625);

  textBox->SetShadowColor(0);
  textBox->SetTextSize(0.025);
}

void DropletFigure::setStyle() {
  setLineStyle();
  setHistStyle();
  setLegendStyle();
}

void DropletFigure::setValues() {
  bulkEdge = dropletPtr->bulk.radius;
  monoEdge = dropletPtr->monolayer.radius;
  dropletHeight = dropletPtr->bulk.height;
  contactAngle = dropletPtr->bulk.contactAngle;
  monoLimits[0] = dropletPtr->monolayer.zlim[0];
  monoLimits[1] = dropletPtr->monolayer.zlim[1];
}

void DropletFigure::setLegendText() {
  // TODO: Do this differently (w/o ss)
  stringstream ss;
  //Add data text box
  ss.str("");
  ss << "Contact angle: " << contactAngle;
  cAText->SetText(0,0,ss.str().data());
  ss.str("");
  ss << "Droplet height: " << dropletHeight;
  dHText->SetText(0,0,ss.str().data());
  ss.str("");
  ss << "Bulk radius: " << bulkEdge;
  bEText->SetText(0,0,ss.str().data());
  ss.str("");
  ss << "Mono radius: " << monoEdge;
  mEText->SetText(0,0,ss.str().data());
}

void DropletFigure::addLegendEntries() {
  legend->AddEntry(gCirclePoints,"Droplet boundary","lp");
  legend->AddEntry(bulkEdgeLine,"Bulk radius","l");
  legend->AddEntry(monoEdgeLine,"Mono radius","l");
  legend->AddEntry(heightLine,"Droplet height","l");
  legend->AddEntry(monoHiLine,"Mono top","l");
  legend->AddEntry(monoLoLine,"Mono bottom","l");
  legend->AddEntry(tangentLine,"Tangent line","l");
}

void DropletFigure::drawHist() {
  hDroplet->Draw("colZ");
}

void DropletFigure::drawLines() {
  // Bulk circle
  circleEllipse = circlePtr->Draw();
  // Contact angle tangent
  tangentLine = circlePtr->DrawTangentLine();
  // Droplet height
  heightLine->SetY1(dropletHeight);
  heightLine->SetY2(dropletHeight);
  heightLine->Draw();

  // Bulk radius
  bulkEdgeLine->SetX1(bulkEdge);
  bulkEdgeLine->SetX2(bulkEdge);
  bulkEdgeLine->Draw();
  // Monolayer radius
  monoEdgeLine->SetX1(monoEdge);
  monoEdgeLine->SetX2(monoEdge);
  monoEdgeLine->Draw();

  // Monolayer
  monoHiLine->SetY1(monoLimits[1]);
  monoHiLine->SetY2(monoLimits[1]);
  monoHiLine->Draw();
  monoLoLine->SetY1(monoLimits[0]);
  monoLoLine->SetY2(monoLimits[0]);
  monoLoLine->Draw();

  //Draw circle points graph
  gCirclePoints->SetMarkerStyle(20);
  gCirclePoints->Draw("same P");
  // Draw circle through graph
  circleEllipse->Draw();
  circleEllipse->SetLineWidth(2);
  circleEllipse->SetFillStyle(0);
}

void DropletFigure::drawLegend() {
  legend->Draw();
  textBox->Draw();
}

void DropletFigure::draw() {
  setValues();
  canvas->cd();
  setStyle();
  drawLines();
  drawLegend();
  canvas->Update();
}

DensFigure::DensFigure(string _title, string _outFile, Droplet &droplet, Substrate &substrate) : Figure(_title, _outFile, *droplet.simDataPtr) {
  dropletPtr = &droplet;
  substratePtr = &substrate;

  hLiquidDens = dropletPtr->hLiquidDens;
  hSubstrateDens = substratePtr->hSubstrateDens;

  // Set bounds
  xlo = simDataPtr->simBounds.zlo;
  xhi = simDataPtr->simBounds.zhi;
  ylo = options.dens_min;
  yhi = options.dens_max;
  cout << "--- yhi = " << yhi << "---" << endl;

  createLines();
  createLegend();
}

DensFigure::~DensFigure() {
  deleteLines();
  deleteLegend();
}

void DensFigure::createLines() {
  monoHiLineDens = new TLine(xlo,ylo,xlo,yhi); //top of monolayer
  monoLoLineDens = new TLine(xlo,ylo,xlo,yhi); //bottom of monolayer
}

void DensFigure::createLegend() {
  legend = new TLegend(.75,.75,.85,.85);
}

void DensFigure::deleteLines() {
  delete monoHiLineDens;
  delete monoLoLineDens;
}

void DensFigure::deleteLegend() {
  delete legend;
}

void DensFigure::setLineStyle() {
  monoHiLineDens->SetLineWidth(3);
  monoHiLineDens->SetLineColor(kRed);
  monoLoLineDens->SetLineWidth(3);
  monoLoLineDens->SetLineColor(kGreen);

  cout << "hLD @ " << hLiquidDens << endl;
  hLiquidDens->SetLineColor(kBlue);
  hLiquidDens->SetLineWidth(2);
  hSubstrateDens->SetLineColor(kOrange+3); //Brown
  hSubstrateDens->SetLineWidth(2);

  // Axis limits
  hLiquidDens->SetAxisRange(xlo, xhi, "X");
  hLiquidDens->SetAxisRange(ylo, yhi, "Y");
  hSubstrateDens->SetAxisRange(xlo, xhi, "X");
  hSubstrateDens->SetAxisRange(ylo, yhi, "Y");

  cout << "Set range to " << ylo << " -> " << yhi << endl;
}

void DensFigure::setLegendStyle() {
  hLiquidDens->SetStats(0);
  hSubstrateDens->SetStats(0);
}

void DensFigure::setStyle() {
  setLineStyle();
  setLegendStyle();
}

void DensFigure::setValues() {
  monoLimits[0] = dropletPtr->monolayer.zlim[0];
  monoLimits[1] = dropletPtr->monolayer.zlim[0];
}

void DensFigure::drawLines() {
  monoHiLineDens->SetX1(monoLimits[1]);
  monoHiLineDens->SetX2(monoLimits[1]);
  // monoHiLineDens->Draw();
  monoLoLineDens->SetX1(monoLimits[0]);
  monoLoLineDens->SetX2(monoLimits[0]);
  // monoLoLineDens->Draw("same");

  cout << "Drawing" << endl;
  hLiquidDens->Draw("hist");
  hSubstrateDens->Draw("hist same"); //Same canvas
  cout << "Drawn" << endl;
}

void DensFigure::drawLegend() {
  legend->AddEntry(hLiquidDens,"Water");
  legend->AddEntry(hSubstrateDens,"Substrate");
  legend->AddEntry(monoLoLineDens,"Mono lower limit","l");
  legend->AddEntry(monoHiLineDens,"Mono upper limit","l");
  legend->Draw();
}

void DensFigure::draw() {
  setValues();
  canvas->cd();
  setStyle();
  drawLines();
  drawLegend();
  canvas->Update();
}

TanhFigure::TanhFigure(string _string, string _outFile, SimData &simData) : Figure(_string, _outFile, simData) {
  // x = z
  xlo = simDataPtr->simBounds.zlo;
  xhi = simDataPtr->simBounds.zhi;
  // y = density
  ylo = simDataPtr->options.dens_max;
  yhi = simDataPtr->options.dens_min;

  createLines();
  createGraph();
  createLegend();
}

TanhFigure::~TanhFigure() {
  deleteLines();
  deleteGraph();
  deleteLegend();
}

void TanhFigure::createLines() {
  //Lines (There are 7)
  ldLine = new TLine(xlo,ylo,xhi,ylo); //Liquid density (horizontal)
  solLine = new TLine(xlo,ylo,xlo,yhi); //Final solution
  x0Line = new TLine(xlo,ylo,xlo,yhi); //x0 final
  x0guessLine = new TLine(xlo,ylo,xlo,yhi); //x0 guess
  lowGuessLine = new TLine(xlo,ylo,xlo,yhi); //Lower edge of boundary
  hiGuessLine = new TLine(xlo,ylo,xlo,yhi); //Upper edge of boundary
  lowBinLine = new TLine(xlo,ylo,xlo,yhi); //Lowest bin used
}

void TanhFigure::createGraph() {
  gPoints = new TGraph();
}

void TanhFigure::createLegend() {
  legend = new TLegend(.65,.65,.85,.85);
  tanhTextBox = new TPaveText();
}

void TanhFigure::deleteGraph() {
  delete gPoints;
}

void TanhFigure::deleteLines() {
  delete ldLine;
  delete solLine;
  delete x0Line;
  delete x0guessLine;
  delete lowGuessLine;
  delete hiGuessLine;
  delete lowBinLine;
}

void TanhFigure::deleteLegend() {
  delete legend;
  delete tanhTextBox;
}

void TanhFigure::addLegendEntries() {
  legend->AddEntry(ldLine,"ld","l");
  legend->AddEntry(solLine,"sol","l");
  legend->AddEntry(x0Line,"x0","l");
  legend->AddEntry(x0guessLine,"x0guess","l");
  legend->AddEntry(lowGuessLine,"lowGuess","l");
  legend->AddEntry(hiGuessLine,"hiGuess","l");
  legend->AddEntry(lowBinLine,"lowBin","l");
}

void TanhFigure::setValues(double _ld, double _w_guess, double _w, double _x0_guess, double _x0) {
  ld = _ld;
  w_guess = _w_guess;
  w = _w;
  x0_guess = _x0_guess;
  x0 = _x0;
}

void TanhFigure::setText() {
  stringstream ss;

  //Text lines (There are 5)
  tanhTexts[0]=tanhTextBox->AddText("ld: ");
  tanhTexts[1]=tanhTextBox->AddText("w_guess: ");
  tanhTexts[2]=tanhTextBox->AddText("w: ");
  tanhTexts[3]=tanhTextBox->AddText("x0_guess: ");
  tanhTexts[4]=tanhTextBox->AddText("x0: ");

  //Text - Where is this bin located?
  // TODO: Delete this?
  // ss.str("");
  // ss << fitType << " pos: " << pos;
  // posText = new TText(.45,.85,ss.str().data());
  // posText->SetNDC();
  // posText->Draw();

  //Set texts
  // TODO: Do without ss
  //pos
  ss.str("");
  ss << "ld: " << ld;
  tanhTexts[0]->SetText(0,0,ss.str().data());
  //w_guess
  ss.str("");
  ss << "w_guess: " << width;
  tanhTexts[1]->SetText(0,0,ss.str().data());
  //w
  ss.str("");
  ss << "w: " << w;
  tanhTexts[2]->SetText(0,0,ss.str().data());
  //x0_guess
  ss.str("");
  ss << "x0_guess: " << x0_guess;
  tanhTexts[3]->SetText(0,0,ss.str().data());
  //x0
  ss.str("");
  ss << "x0: " << x0;
  tanhTexts[4]->SetText(0,0,ss.str().data());
  ss.str("");
}

void TanhFigure::setLinePosition() {
  //Set lines to span plot box
  //ld (horizontal)
  tanhLines[0]->SetY1(ld);
  tanhLines[0]->SetY2(ld);
  tanhLines[0]->Draw();
  //x0
  tanhLines[2]->SetX1(x0);
  tanhLines[2]->SetX2(x0);
  tanhLines[3]->Draw();
  //x0guess
  tanhLines[3]->SetX1(x0_guess);
  tanhLines[3]->SetX2(x0_guess);
  tanhLines[3]->Draw();
  //lowGuess
  tanhLines[4]->SetX1(lowGuess);
  tanhLines[4]->SetX2(lowGuess);
  tanhLines[4]->Draw();
  //hiGuess
  tanhLines[5]->SetX1(hiGuess);
  tanhLines[5]->SetX2(hiGuess);
  tanhLines[5]->Draw();
  //lowBin
  tanhLines[6]->SetX1(startPoint);
  tanhLines[6]->SetX2(startPoint);
  tanhLines[6]->Draw();
  //sol
  tanhLines[1]->SetX1(val);
  tanhLines[1]->SetX2(val);
  tanhLines[1]->Draw();

  // TODO: Use variable instead of 7
  for(int i=1;i<7;i++)
    {
      //Vertical lines have identical x coordinates
      if(tanhLines[i]->GetX1()==tanhLines[i]->GetX2())
        {
          tanhLines[i]->SetY1(ylo);
          tanhLines[i]->SetY2(yhi);
        }
      //Horizontal lines do not
      else
        {
          tanhLines[i]->SetX1(xlo);
          tanhLines[i]->SetX2(xhi);
        }
    }
}

void TanhFigure::setLineStyle() {
  ldLine->SetLineColor(kYellow);
  ldLine->SetLineWidth(3);
  solLine->SetLineColor(kBlack);
  solLine->SetLineWidth(3);
  x0Line->SetLineColor(kBlue);
  x0Line->SetLineWidth(3);
  x0guessLine->SetLineColor(kCyan);
  x0guessLine->SetLineWidth(3);
  lowGuessLine->SetLineColor(kViolet);
  lowGuessLine->SetLineWidth(3);
  hiGuessLine->SetLineColor(kGreen);
  hiGuessLine->SetLineWidth(3);
  lowBinLine->SetLineColor(kOrange+7);
  lowBinLine->SetLineWidth(3);
}

void TanhFigure::setPointStyle() {
  gPoints->SetMarkerStyle(20);
}

void TanhFigure::setLegendStyle() {
  //Text to show data
  tanhTextBox->SetX1NDC(.65);
  tanhTextBox->SetY1NDC(.45);
  tanhTextBox->SetX2NDC(.85);
  tanhTextBox->SetY2NDC(.625);

  tanhTextBox->SetShadowColor(0);
  tanhTextBox->SetTextSize(0.025);
}

void TanhFigure::setStyle() {
  setLineStyle();
  setPointStyle();
  setLegendStyle();
}

void TanhFigure::drawLines() {
  for(int i=0; i<7; i++) {
    tanhLines[i]->Draw();
  }

  //Draw solution
  TLine *solLine= new TLine(val,0,val,1);
  solLine->SetLineWidth(3);
  solLine->SetLineColor(kOrange);
  solLine->Draw();
}

void TanhFigure::drawPoints() {
  gPoints->Draw("same p");
}

void TanhFigure::drawLegend() {
  legend->Draw();
}

void TanhFigure::draw() {
  canvas->cd();
  setStyle();
  drawLines();
  drawPoints();
  drawLegend();
  canvas->Update();
}

// TODO: What is this?
TGraph *horizontalHist(TH1D* hist)
{
    char *title, *xLabel, *yLabel;
    int n;
    double *x,*y, xlo, xhi;

    //Copy points
    TGraph* g = new TGraph(hist);
    n=g->GetN();
    x=g->GetX();
    y=g->GetY();
    xlo=hist->GetXaxis()->GetXmin();
    xhi=hist->GetXaxis()->GetXmax();

    //Copy title and axis labels
    title=(char*) hist->GetTitle();
    xLabel=(char*) hist->GetXaxis()->GetTitle();
    yLabel=(char*) hist->GetYaxis()->GetTitle();

    //Flip x & y
    TGraph* g1 = new TGraph(n,y,x);

    //Copy titles
    g1->SetTitle(title);
    g1->GetXaxis()->SetTitle(yLabel);
    g1->GetYaxis()->SetTitle(xLabel);

    //Copy line attributes
    g1->SetLineStyle(g->GetLineStyle());
    g1->SetLineWidth(g->GetLineWidth());
    g1->SetLineColor(g->GetLineColor());

    //Copy fill attributes
    g1->SetFillStyle(g->GetFillStyle());
    g1->SetFillColor(g->GetFillColor());

    //Draw
    g1->Draw("AL");
    g1->GetYaxis()->SetRangeUser(xlo,xhi);
    g1->Draw("AL");

    //Delete
    delete g;
    return g1;
}
