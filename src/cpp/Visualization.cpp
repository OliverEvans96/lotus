
#include "Visualization.h"

///////////////////
// Visualization //
///////////////////

Figure::Figure(string _title, string _outFile, SimData* _simDataPtr) {
  title = _title;
  outFile = _outFile;
  options = _options;
  simDataPtr = _simDataPtr;
  boundsPtr = simDataPtr->simBounds;

  createCanvas(name);
  setCanvasStyle();
}

Figure::~Figure() {
  delete canvas;
}

void Figure::createCanvas() {
  width = options.plot_width;
  height = (int) round(width / options.plot_aspect);
  canvas = new TCanvas(name, name, width, height);
}

void Figure::setCanvasStyle() {
  // Use OpenGL for antialiasing
  gStyle->SetCanvasPreferGL(true);
}

void Figure::saveImage() {
  canvas->SaveAs(outFile);
}

void Figure::saveROOT() {
  // TODO: save w/ .C extension
}

void Figure::save() {
  if(options.saveImages) {
    saveImage();
  }

  if(options.saveROOT) {
    saveROOT();
  }
}

// TODO: Constructor
DropletFigure::DropletFigure(TH2D* _hDroplet, string _title, string _outFile, SimData* _simDataPtr) : Figure(_title, _outFile, _simDataPtr) {
  hDroplet = _hDroplet;
  xlo = 0;
  xhi = simDataPtr->plot_rmax;
  ylo = 0;
  yhi = simDataPtr->plot_zmax;
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
  hist->SetMinimum(colzMin);
  hist->SetMaximum(colzMax);
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

void DropletFigure::drawHist() {
  hist->Draw("colZ")
}

void DropletFigure::drawLines() {
  // Bulk circle
  circleEllipse = circle.Draw();
  // Contact angle tangent
  tangentLine = circle.DrawTangentLine();
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
  circlePointsGraph->SetMarkerStyle(20);
  circlePointsGraph->Draw("same P");
}

void DropletFigure::drawLegend() {
  hALegend->Draw();
  textBox->Draw();
}

void DropletFigure::setLegendText() {
  // TODO: Do this differently (w/o ss)
  //Add data text box
  title.str("");
  title << "Contact angle: " << contactAngle;
  cAText->SetText(0,0,title.str().data());
  title.str("");
  title << "Droplet height: " << dropletHeight;
  dHText->SetText(0,0,title.str().data());
  title.str("");
  title << "Bulk radius: " << bulkEdge;
  bEText->SetText(0,0,title.str().data());
  title.str("");
  title << "Mono radius: " << monoEdge;
  mEText->SetText(0,0,title.str().data());
}

void DropletFigure::addLegendEntries() {
  hALegend->AddEntry(circlePointsGraph,"Droplet boundary","lp");
  hALegend->AddEntry(bulkEdgeLine,"Bulk radius","l");
  hALegend->AddEntry(monoEdgeLine,"Mono radius","l");
  hALegend->AddEntry(heightLine,"Droplet height","l");
  hALegend->AddEntry(monoHiLine,"Mono top","l");
  hALegend->AddEntry(monoLoLine,"Mono bottom","l");
  hALegend->AddEntry(tangentLine,"Tangent line","l");
}

void DropletFigure::drawLegend() {
  //2D Density Hist Legend
}

void DropletFigure::draw() {
  drawLines();
  drawLegend();
}

DensFigure::DensFigure() {
  createLines();
  createLegend();
}

DensFigure::~DensFigure() {
  deleteLines();
  deleteLegend();
}

void DensFigure::createLines() {
}

void DensFigure::createLegend() {
  densLeg = new TLegend(.75,.75,.85,.85);
}

void DensFigure::deleteLines() {
}

void DensFigure::deleteLegend() {
  delete densLeg;
}

void DensFigure::setLineStyle() {
  monoHiLineDens->SetLineWidth(3);
  monoHiLineDens->SetLineColor(kRed);
  monoLoLineDens->SetLineWidth(3);
  monoLoLineDens->SetLineColor(kGreen);

  hLiquidDens->SetLineColor(kBlue);
  hLiquidDens->SetLineWidth(2);
  hSubstrateDens->SetLineColor(kOrange+3); //Brown
  hSubstrateDens->SetLineWidth(2);

  // TODO: Set limits from options
  // Axis limits
  hLiquidDens->GetYaxis()->SetRangeUser(0,6);
}


void DensFigure::setLegendStyle() {}

void DensFigure::setStyle() {}

void DensFigure::drawLines() {
  monoHiLineDens->SetX1(monoLimits[1]);
  monoHiLineDens->SetX2(monoLimits[1]);
  monoHiLineDens->Draw();
  monoLoLineDens->SetX1(monoLimits[0]);
  monoLoLineDens->SetX2(monoLimits[0]);
  monoLoLineDens->Draw();

  hLiquidDens->Draw();
  hSubstrateDens->Draw("SAME"); //Same canvas
}
void DensFigure::drawLegend() {
  densLeg->AddEntry(hLiquidDens,"Water");
  densLeg->AddEntry(hSubstrateDens,"Substrate");
  densLeg->AddEntry(monoLoLineDens,"Mono lower limit","l");
  densLeg->AddEntry(monoHiLineDens,"Mono upper limit","l");
  densLeg->Draw();
}

// TODO: Remove
void DensFigure::draw() {
  canvas->cd();
  drawLines();
  drawLegend();
}

// TODO: Modify for DensFigure from Droplet
void DensFigure::plotDensity(char* filename) {
  stringstream ss;
  ss << options.outLoc << "/" << filename;
  cDroplet->cd();
  hDroplet->Draw("colZ");
  if(options.verbose) {
    cout << "Saving droplet hist @ '" << ss.str().data() << "'" << endl;
  }
}

TanhFigure::TanhFigure() {
  // TODO: set xlo, xhi, etc.
  // x = z
  xlo = simDataPtr->simBounds.zlo;
  xhi = simDataPtr->simBounds.zhi;
  // y = density
  ylo = simDataPtr->options.dens_max;
  yhi = simDataPtr->options.dens_min;

  createLines();
}

TanhFigure::~TanhFigure() {
  delete monoHiLineDens;
  delete monoLoLineDens;
}

void TanhFigure::createLines() {
  monoHiLineDens = new TLine(xlo,ylo,xlo,yhi); //top of monolayer
  monoLoLineDens = new TLine(xlo,ylo,xlo,yhi); //bottom of monolayer
}

void TanhFigure::setStyle() {
}

void TanhFigure::fillLegend() {}
void TanhFigure::fillPaveText() {}
void TanhFigure::draw() {}


/////////////

// from Droplet::findBoundaryPoints

// TCanvas* c5 = new TCanvas();
//  c5->cd();
// circlePointsGraph->Draw("APL");
//  c5->SaveAs("test.C");
//  delete c5;

///////////////

// from CircularBulk::fitCircle

//TF2 *circleFitFunction = new TF2("circleFitFunction","(x-[0])*(x-[0])+(y-[1])*(y-[1])-[2]*[2]",0,120,25,100);

//TF1 *circleFitFunction = new TF1("circleFitFunction","sqrt([0]*[0]-(x-[1])*(x-[1]))+[2]",0,120);
//circleFitFunction->GetXaxis()->SetRangeUser(0,xMax);
//circleFitFunction->GetYaxis()->SetRangeUser(30,60);
//circleFitFunction->SetParameters(-13,-134,2200);

//cout << endl << "GRAPH JUST BEFORE FIT" << endl;
//for(int i=0;i<gCirclePoints->GetN();i++)
//    cout << gCirclePoints->GetX()[i] << " " << gCirclePoints->GetY()[i] << endl; //" " << gCirclePoints->GetZ()[i] << endl;

///

//gCirclePoints->Fit("circleFitFunction","RW");

//double x0=circleFitFunction->GetParameter(0);
//double y0=circleFitFunction->GetParameter(1);
//double r=circleFitFunction->GetParameter(2);
//TEllipse *e = new TEllipse(x0,y0,r,r);
//e->SetLineWidth(2);
//e->SetFillStyle(0);
//e->Draw();

///

//TCanvas *cCirc = new TCanvas();
//cCirc->cd();
//cout << "CURRENT CANVAS: " << gPad->GetCanvas()->GetName() << endl;
//gCirclePoints->Draw("AL");
//e->Draw();
//cCirc->SaveAs("circleFit.C");
//stringstream s;
//s << "img/circles/step" << timestep << ".png";
//cCirc->SaveAs(s.str().data());

//////////////////////

// from CircularBulk::solveTanhFit

    //if(frameStep==8810000)
    //draw=true;
    //else
    //    draw=false;

    // //Save image
    // if(draw && val>0)
    // {
    //     //Canvas
    //     TCanvas *cTanh = new TCanvas();

    //     //Title
    //     stringstream title;

    //     //Draw hist
    //     hist->GetXaxis()->SetRangeUser(tanhRangeUser[0],tanhRangeUser[1]);
    //     hist->GetYaxis()->SetRangeUser(tanhRangeUser[2],tanhRangeUser[3]);
    //     hist->Draw();

    //     //Text - Where is this bin located?
    //     title.str("");
    //     title << fitType << " pos: " << pos;
    //     TText *posText = new TText(.45,.85,title.str().data());
    //     posText->SetNDC();
    //     posText->Draw();

    //     //Set annotation lines

    //     //Set lines to span plot box
    //     for(int i=1;i<7;i++)
    //     {
    //         //Vertical lines have identical x coordinates
    //         if(tanhLines[i]->GetX1()==tanhLines[i]->GetX2())
    //         {
    //             tanhLines[i]->SetY1(tanhRangeUser[2]);
    //             tanhLines[i]->SetY2(tanhRangeUser[3]);
    //         }
    //         //Horizontal lines do not
    //         else
    //         {
    //             tanhLines[i]->SetX1(tanhRangeUser[0]);
    //             tanhLines[i]->SetX2(tanhRangeUser[1]);
    //         }
    //     }

    //     //ld (horizontal)
    //     tanhLines[0]->SetY1(ld);
    //     tanhLines[0]->SetY2(ld);
    //     tanhLines[0]->Draw();
    //     //x0
    //     tanhLines[2]->SetX1(c);
    //     tanhLines[2]->SetX2(c);
    //     tanhLines[3]->Draw();
    //     //x0guess
    //     tanhLines[3]->SetX1(boundary);
    //     tanhLines[3]->SetX2(boundary);
    //     tanhLines[3]->Draw();
    //     //lowGuess
    //     tanhLines[4]->SetX1(lowGuess);
    //     tanhLines[4]->SetX2(lowGuess);
    //     tanhLines[4]->Draw();
    //     //hiGuess
    //     tanhLines[5]->SetX1(hiGuess);
    //     tanhLines[5]->SetX2(hiGuess);
    //     tanhLines[5]->Draw();
    //     //lowBin
    //     tanhLines[6]->SetX1(startPoint);
    //     tanhLines[6]->SetX2(startPoint);
    //     tanhLines[6]->Draw();
    //     //sol
    //     tanhLines[1]->SetX1(val);
    //     tanhLines[1]->SetX2(val);
    //     tanhLines[1]->Draw();

    //     //Set texts
    //     //pos
    //     title.str("");
    //     title << "pos: " << pos;
    //     tanhTexts[0]->SetText(0,0,title.str().data());
    //     //w_guess
    //     title.str("");
    //     title << "w_guess: " << width;
    //     tanhTexts[1]->SetText(0,0,title.str().data());
    //     //w
    //     title.str("");
    //     title << "w: " << w;
    //     tanhTexts[2]->SetText(0,0,title.str().data());
    //     //x0_guess
    //     title.str("");
    //     title << "x0_guess: " << boundary;
    //     tanhTexts[3]->SetText(0,0,title.str().data());
    //     //x0
    //     title.str("");
    //     title << "x0: " << c;
    //     tanhTexts[4]->SetText(0,0,title.str().data());
    //     title.str("");

    //     //Draw annotations
    //     tanhTextBox->Draw();
    //     tanhLegend->Draw();

    //     //Draw graph
    //     tanhPointsGraph->SetMarkerStyle(20);
    //     tanhPointsGraph->Draw("same p");

    //     /*
    //     //Draw solution
    //     TLine *solLine= new TLine(val,0,val,1);
    //     solLine->SetLineWidth(3);
    //     solLine->SetLineColor(kOrange);
    //     solLine->Draw();
    //     */

    //     //Title
    //     title.str("");
    //     title << "img/tanh/step" << frameStep << "_" << fitType << setw(3) << setfill('0') << fitNum << ".png";
    //     hist->SetTitle(title.str().data());

    //     //Save
    //     cTanh->SaveAs(title.str().data());

    //     delete cTanh;
    // }

//////////////////////////////////



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
