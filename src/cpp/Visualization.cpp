
#include "Visualization.h"

///////////////////
// Visualization //
///////////////////

Figure::Figure(string _title, SimData &simData) {
  title = _title;
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

// TODO: Save relative to image directory
void Figure::save(char* filename) {
  canvas->SaveAs(filename);
}

DropletFigure::DropletFigure(string _title, Droplet &droplet) : Figure(_title, *droplet.simDataPtr){
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
  cout << "setting ML0 = " << dropletPtr->monolayer.zlim[0] << endl;
  monoLimits[0] = dropletPtr->monolayer.zlim[0];
  cout << "setting ML1 = " << dropletPtr->monolayer.zlim[1] << endl;
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

DensFigure::DensFigure(string _title, Droplet &droplet, Substrate &substrate) : Figure(_title, *droplet.simDataPtr) {
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

void DensFigure::drawHists() {
  cout << "Drawing" << endl;
  hLiquidDens->Draw("hist same");
  hSubstrateDens->Draw("hist same"); //Same canvas
  cout << "Drawn" << endl;
}

void DensFigure::drawLines() {
  monoHiLineDens->SetX1(monoLimits[1]);
  monoHiLineDens->SetX2(monoLimits[1]);
  monoHiLineDens->Draw("same");
  monoLoLineDens->SetX1(monoLimits[0]);
  monoLoLineDens->SetX2(monoLimits[0]);
  monoLoLineDens->Draw("same");
}

void DensFigure::drawLegend() {
  legend->AddEntry(hLiquidDens,"Liquid");
  legend->AddEntry(hSubstrateDens,"Substrate");
  legend->AddEntry(monoLoLineDens,"Mono lower limit","l");
  legend->AddEntry(monoHiLineDens,"Mono upper limit","l");
  legend->Draw("same");
}

void DensFigure::draw() {
  setValues();
  canvas->cd();
  setStyle();
  drawHists();
  drawLines();
  drawLegend();
  canvas->Update();
}

TanhFigure::TanhFigure(string _title, TanhFit &tanhFit) : Figure(_title, *tanhFit.simDataPtr) {
  tanhFitPtr = &tanhFit;

  // x = z
  xlo = simDataPtr->simBounds.zlo;
  xhi = simDataPtr->simBounds.zhi;
  // y = density
  ylo = simDataPtr->options.dens_max;
  yhi = simDataPtr->options.dens_min;

  createLines();
  // createGraph();
  createLegend();
  addLegendEntries();
}

TanhFigure::~TanhFigure() {
  deleteLines();
  // deleteGraph();
  deleteLegend();
}

void TanhFigure::createLines() {
  //Lines (There are 7)
  ldLine = new TLine(xlo,ylo,xhi,ylo); // Liquid density (horizontal)
  halfLdLine = new TLine(xlo,ylo,xhi,ylo); // HalfLiquid density (horizontal)
  x0Line = new TLine(xlo,ylo,xlo,yhi); // x0 final

  tanhLines[0] = ldLine;
  tanhLines[1] = halfLdLine;
  tanhLines[2] = x0Line;
}

void TanhFigure::createGraph() {
  gPoints = new TGraph();
}

void TanhFigure::createLegend() {
  legend = new TLegend(.65,.65,.85,.85);
  tanhTextBox = new TPaveText();
}

void TanhFigure::deleteLines() {
  delete ldLine;
  delete halfLdLine;
  delete x0Line;
}

void TanhFigure::deleteGraph() {
  delete gPoints;
}

void TanhFigure::deleteLegend() {
  delete legend;
  delete tanhTextBox;
}

void TanhFigure::addLegendEntries() {
  legend->AddEntry(ldLine,"ld","l");
  legend->AddEntry(halfLdLine,"ld/2","l");
  legend->AddEntry(x0Line,"x0","l");
}

void TanhFigure::setValues() {
  ld = tanhFitPtr->getLiquidDensity();
  w = tanhFitPtr->getWidth();
  x0 = tanhFitPtr->getBoundary();
}

void TanhFigure::setText() {
  stringstream ss;

  //Text lines (There are 5)
  tanhTexts[0]=tanhTextBox->AddText("ld: ");
  tanhTexts[2]=tanhTextBox->AddText("w: ");
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
  //w
  ss.str("");
  ss << "w: " << w;
  tanhTexts[2]->SetText(0,0,ss.str().data());
  //x0
  ss.str("");
  ss << "x0: " << x0;
  tanhTexts[4]->SetText(0,0,ss.str().data());
  ss.str("");
}

void TanhFigure::setLinePositions() {
  //Set lines to span plot box
  //ld (horizontal)
  ldLine->SetY1(ld);
  ldLine->SetY2(ld);
  //ld/2 (horizontal)
  halfLdLine->SetY1(ld/2);
  halfLdLine->SetY2(ld/2);
  //x0
  x0Line->SetX1(x0);
  x0Line->SetX2(x0);

  for(int i=1;i<numLines;i++)
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

void TanhFigure::setHistStyle() {
  tanhFitPtr->hTanh->SetStats(0);
}

void TanhFigure::setLineStyle() {
  ldLine->SetLineColor(kYellow);
  ldLine->SetLineWidth(3);
  halfLdLine->SetLineColor(kBlack);
  halfLdLine->SetLineWidth(3);
  x0Line->SetLineColor(kBlue);
  x0Line->SetLineWidth(3);
}

void TanhFigure::setGraphStyle() {
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
  setHistStyle();
  setLineStyle();
  // setGraphStyle();
  setLegendStyle();
}

void TanhFigure::drawHist() {
  tanhFitPtr->hTanh->Draw("hist");

  for(int i=1; i<=tanhFitPtr->hTanh->GetXaxis()->GetNbins(); i++) {
    cout << "hTanh " << i << " " << tanhFitPtr->hTanh->GetBinContent(i) << endl;
  }
}

void TanhFigure::drawFunction() {
  tanhFitPtr->fTanh->Draw("same");
}

void TanhFigure::drawLines() {
  for(int i=0; i<numLines; i++) {
    tanhLines[i]->Draw("same");
  }
}

void TanhFigure::drawGraph() {
  gPoints->Draw("same p");
}

void TanhFigure::drawLegend() {
  legend->Draw("same");
}

void TanhFigure::draw() {
  setValues();
  setText();
  setLinePositions();

  canvas->cd();
  setStyle();

  drawHist();
  drawFunction();
  drawLines();
  // drawGraph();
  drawLegend();
  canvas->Update();
}
