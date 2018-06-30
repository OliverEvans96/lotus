
#include "Visualization.h"

///////////////////
// Visualization //
///////////////////

// Create figure. title must be unique among all Figures.
Figure::Figure(const string _title, SimData &simData) {
  title = _title;
  simDataPtr = &simData;
  options = simDataPtr->options;

  setOutputDir(options.outLoc.data());
  createDirs();
  createCanvas();
  setCanvasStyle();
}

Figure::~Figure() {
  delete canvas;
}

void Figure::createCanvas() {
  char str[256];
  width = options.plot_width;
  height = (int) round(width / options.plot_aspect);
  if(options.verbose) {
    cout << "canvas width = " << width << endl;
    cout << "canvas height = " << height << endl;
  }
  strcpy(str, title.data());
  canvas = new TCanvas(str, str, width, height);
}

void Figure::setCanvasStyle() {
  // Use OpenGL for antialiasing
  gStyle->SetCanvasPreferGL(true);
}

void Figure::setOutputDir(const char* path) {
  strcpy(outputDir, path);
}

// Get filename, composed from title and timestep
void Figure::getFilename(char* filename, const char* suffix) {
  // suffix determines file type (e.g. "png", "C")
  sprintf(filename, "%09d.%s", simDataPtr->framePtr->time, suffix);
}

// Create `outputDir`/`subdir`/`title`
void Figure::createDirBase(const char* subdir) {
  char outerDir[256];
  char innerDir[256];
  // `outputDir`
  if(!dir_exists(outputDir)) {
    mkdir(outputDir, S_IRWXU);
  }

 // `outputDir`/`subdir`
  joinPath(outerDir, outputDir, subdir);
  if(!dir_exists(outerDir)) {
    mkdir(outerDir, S_IRWXU);
  }

  // `outputDir`/`subdir`/`title`
  joinPath(innerDir, outerDir, title.data());
  if(!dir_exists(innerDir)) {
    mkdir(innerDir, S_IRWXU);
  }
}

void Figure::createImageDir() {
  createDirBase("img");
}

void Figure::createROOTDir() {
  createDirBase("root");
}

void Figure::createDirs() {
  if(options.saveImages)
    createImageDir();
  if(options.saveROOT)
    createROOTDir();
}

// Save to `outputDir`/`subdir`/`title`/`timestep`.`suffix`
void Figure::saveBase(const char* subdir, const char* suffix) {
  char filename[256];
  char outerDir[256];
  char innerDir[256];
  char path[256];

  getFilename(filename, suffix);
  joinPath(outerDir, outputDir, subdir);
  joinPath(innerDir, outerDir, title.data());
  joinPath(path, innerDir, filename);
  canvas->SaveAs(path);
  if(options.verbose)
    cout << "saved canvas " << canvas << " @ '" << path << "'" << endl;
}


void Figure::saveImage() {
  saveBase("img", "png");
}

void Figure::saveROOT() {
  saveBase("root", "C");
}

void Figure::save() {
  if(options.saveImages)
    saveImage();
  if(options.saveROOT)
    saveROOT();
}

DropletFigure::DropletFigure(const string _title, Droplet &droplet) : Figure(_title, *droplet.simDataPtr){
  dropletPtr = &droplet;
  hDroplet = dropletPtr->hDroplet;
  gCirclePoints = dropletPtr->bulk.gCirclePoints;
  circlePtr = &dropletPtr->bulk.circle;

  xlo = 0;
  xhi = options.plot_rmax;
  ylo = 0;
  yhi = options.plot_zmax;

  createLines();
  createCircle();
  createLegend();
}

DropletFigure::~DropletFigure() {
  deleteLines();
  deleteCircle();
  deleteLegend();
}

void DropletFigure::createLines() {
  if(options.fitCircle) {
    bulkEdgeLine = new TLine(xlo,ylo,xlo,yhi); //Bulk edge
    heightLine = new TLine(xlo,ylo,xhi,ylo); //Droplet height
    tangentLine = new TLine(xlo,ylo,xhi,ylo); //Contact angle
  }
  if(options.monolayer) {
    monoEdgeLine = new TLine(xlo,ylo,xlo,yhi); //Mono edge
    monoHiLine = new TLine(xlo,ylo,xhi,ylo); //top of monolayer
    monoLoLine = new TLine(xlo,ylo,xhi,ylo); //bottom of monolayer
  }
}

void DropletFigure::createCircle() {
  if(options.fitCircle) {
    eCircle = new TEllipse();
    eGuessCircle = new TEllipse();
  }
}

void DropletFigure::createLegend() {
  legend = new TLegend(.65,.65,.85,.85);
  textBox = new TPaveText();
  addLegendEntries();
}

void DropletFigure::deleteLines() {
  if(options.fitCircle) {
    delete bulkEdgeLine;
    delete heightLine;
    delete tangentLine;
  }

  if(options.monolayer) {
    delete monoEdgeLine;
    delete monoHiLine;
    delete monoLoLine;
  }
}

void DropletFigure::deleteCircle() {
  if(options.fitCircle) {
    delete eCircle;
    delete eGuessCircle;
  }
}

void DropletFigure::deleteLegend() {
  delete legend;
  delete textBox;
}

void DropletFigure::setTitle() {
  char title[256];
  // Note timestep in title
  sprintf(title, "%s, t=%09d", simDataPtr->options.dumpfile.data(), simDataPtr->framePtr->time);
  hDroplet->SetTitle(title);
}

void DropletFigure::setAxisLabels() {
  hDroplet->SetXTitle("r (#AA)");
  hDroplet->SetYTitle("z (#AA)");
  hDroplet->GetXaxis()->CenterTitle();
  hDroplet->GetYaxis()->CenterTitle();
}

void DropletFigure::setLineStyle() {
  //Line properties
  if(options.fitCircle) {
    bulkEdgeLine->SetLineWidth(3);
    bulkEdgeLine->SetLineColor(kGreen);
    heightLine->SetLineWidth(3);
    heightLine->SetLineColor(kOrange+3); // Brown
    tangentLine->SetLineColor(kViolet);
    tangentLine->SetLineWidth(3);
  }

  if(options.monolayer) {
    monoEdgeLine->SetLineWidth(3);
    monoEdgeLine->SetLineColor(kRed);

    monoHiLine->SetLineWidth(3);
    monoHiLine->SetLineColor(kBlue);
    monoHiLine->SetLineStyle(9);
    monoLoLine->SetLineWidth(3);
    monoLoLine->SetLineColor(kBlue);
    monoLoLine->SetLineStyle(2);
  }
}

void DropletFigure::setHistStyle() {
  hDroplet->SetMinimum(0.0);
  hDroplet->SetMaximum(options.densMax);
}

void DropletFigure::setCircleStyle() {
  if(options.fitCircle) {
    eCircle->SetLineWidth(2);
    eCircle->SetFillStyle(0);

    eGuessCircle->SetLineWidth(2);
    eGuessCircle->SetFillStyle(0);
    eGuessCircle->SetLineColor(kBlue);
  }
}

void DropletFigure::setGraphStyle() {
  gCirclePoints->SetMarkerStyle(20);
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
  setHistStyle();
  setLineStyle();
  setCircleStyle();
  setGraphStyle();
  setLegendStyle();
}

void DropletFigure::setValues() {
  if(options.fitCircle) {
    bulkEdge = dropletPtr->bulk.radius;
    monoEdge = dropletPtr->monolayer.radius;
    dropletHeight = dropletPtr->bulk.height;
    contactAngle = dropletPtr->bulk.contactAngle;
  }

  if(options.monolayer) {
    monoLimits[0] = dropletPtr->monolayer.zlim[0] - simDataPtr->substrateTop;
    monoLimits[1] = dropletPtr->monolayer.zlim[1] - simDataPtr->substrateTop;
  }
}

void DropletFigure::setLegendText() {
  char str[256];
  if(options.fitCircle) {
    sprintf(str, "Contact angle: %6.2f#circ", contactAngle);
    cAText->SetText(0,0,str);
    sprintf(str, "Droplet height: %6.2f", dropletHeight);
    dHText->SetText(0,0,str);
    sprintf(str, "Bulk radius: %6.2f", bulkEdge);
    bEText->SetText(0,0,str);
  }

  if(options.monolayer) {
    sprintf(str, "Mono radius: %6.2f", monoEdge);
    mEText->SetText(0,0,str);
    sprintf(str, "Mono bottom: %6.2f", monoLimits[0]);
    mLText->SetText(0,0,str);
    sprintf(str, "Mono top: %6.2f", monoLimits[1]);
    mHText->SetText(0,0,str);
  }
}

void DropletFigure::addLegendEntries() {
  legend->AddEntry(gCirclePoints,"Droplet boundary","lp");

  if(options.fitCircle) {
    legend->AddEntry(bulkEdgeLine, "bulk radius", "l");
    legend->AddEntry(heightLine, "height", "l");
    legend->AddEntry(tangentLine,"Tangent line","l");

    cAText = textBox->AddText("Contact angle");
    dHText = textBox->AddText("Droplet height");
    bEText = textBox->AddText("Bulk radius");
  }

  if(options.monolayer) {
    legend->AddEntry(monoEdgeLine, "mono radius", "l");
    legend->AddEntry(monoHiLine, "mono top", "l");
    legend->AddEntry(monoLoLine, "mono bottom", "l");

    mEText = textBox->AddText("Mono radius");
    mHText = textBox->AddText("Mono top");
    mLText = textBox->AddText("mLText");
  }
}

void DropletFigure::drawHist() {
  hDroplet->Draw("colZ");
}

void DropletFigure::drawLines() {
  if(options.fitCircle) {
    circlePtr->SetTangentLine(tangentLine);
    tangentLine->Draw("same");

    // Droplet height
    heightLine->SetY1(dropletHeight);
    heightLine->SetY2(dropletHeight);
    heightLine->Draw("same");

    // Bulk radius
    bulkEdgeLine->SetX1(bulkEdge);
    bulkEdgeLine->SetX2(bulkEdge);
    bulkEdgeLine->Draw("same");
  }

  if(options.monolayer) {
    // Monolayer radius
    monoEdgeLine->SetX1(monoEdge);
    monoEdgeLine->SetX2(monoEdge);
    monoEdgeLine->Draw("same");

    // Monolayer
    monoHiLine->SetY1(monoLimits[1]);
    monoHiLine->SetY2(monoLimits[1]);
    monoHiLine->Draw("same");
    monoLoLine->SetY1(monoLimits[0]);
    monoLoLine->SetY2(monoLimits[0]);
    monoLoLine->Draw("same");
  }
}

void DropletFigure::drawCircle() {
  eCircle->SetR1(circlePtr->getRadius());
  eCircle->SetR2(circlePtr->getRadius());
  eCircle->SetX1(circlePtr->getXCenter());
  eCircle->SetY1(circlePtr->getYCenter());
  eCircle->Draw("same");

  eGuessCircle->SetR1(circlePtr->gr);
  eGuessCircle->SetR2(circlePtr->gr);
  eGuessCircle->SetX1(circlePtr->gx0);
  eGuessCircle->SetY1(circlePtr->gy0);
  // eGuessCircle->Draw("same");
}

void DropletFigure::drawGraph() {
  gCirclePoints->Draw("same P");
}

void DropletFigure::drawLegend() {
  legend->Draw("same");
  textBox->Draw("same");
}

void DropletFigure::draw() {
  setValues();
  canvas->cd();
  setStyle();
  setLegendText();
  setTitle();
  setAxisLabels();
  drawHist();
  drawLines();
  if(options.fitCircle) {
    drawCircle();
  }
  drawGraph();
  drawLegend();
  canvas->Update();
}

DensFigure::DensFigure(const string _title, Droplet &droplet, Substrate &substrate) : Figure(_title, *droplet.simDataPtr) {
  dropletPtr = &droplet;
  hLiquidDens = dropletPtr->hLiquidDens;

  if(options.substrate) {
    substratePtr = &substrate;
    hSubstrateDens = substratePtr->hSubstrateDens;
  }

  // Set bounds
  xlo = simDataPtr->simBounds.zlo;
  xhi = simDataPtr->simBounds.zhi;
  ylo = options.dens_min;
  yhi = options.dens_max;

  createLines();
  createLegend();
}

DensFigure::~DensFigure() {
  deleteLines();
  deleteLegend();
}

void DensFigure::createLines() {
  if(options.monolayer) {
    monoHiLineDens = new TLine(xlo,ylo,xlo,yhi); //top of monolayer
    monoLoLineDens = new TLine(xlo,ylo,xlo,yhi); //bottom of monolayer
  }
}

void DensFigure::createLegend() {
  legend = new TLegend(.65,.65,.85,.85);

  legend->AddEntry(hLiquidDens,"Liquid");

  if(options.substrate) {
    legend->AddEntry(hSubstrateDens,"Substrate");
  }

  if(options.monolayer) {
    legend->AddEntry(monoLoLineDens,"Mono bottom","l");
    legend->AddEntry(monoHiLineDens,"Mono top","l");
  }
}

void DensFigure::deleteLines() {
  if(options.monolayer) {
    delete monoHiLineDens;
    delete monoLoLineDens;
  }
}

void DensFigure::deleteLegend() {
  delete legend;
}

void DensFigure::setTitle() {
  char title[256];
  // Note timestep in title
  sprintf(title, "%s, t=%09d", simDataPtr->options.dumpfile.data(), simDataPtr->framePtr->time);
  hLiquidDens->SetTitle(title);
}

void DensFigure::setAxisLabels() {
  hLiquidDens->SetXTitle("z (#AA)");
  hLiquidDens->SetYTitle("mass density (g/cc)");
  hLiquidDens->GetXaxis()->CenterTitle();
  hLiquidDens->GetYaxis()->CenterTitle();
}

void DensFigure::setLineStyle() {
  if(options.monolayer) {
    monoHiLineDens->SetLineWidth(3);
    monoHiLineDens->SetLineColor(kRed);
    monoLoLineDens->SetLineWidth(3);
    monoLoLineDens->SetLineColor(kGreen);
  }

  hLiquidDens->SetLineColor(kBlue);
  hLiquidDens->SetLineWidth(2);
  hLiquidDens->SetAxisRange(xlo, xhi, "X");
  hLiquidDens->SetAxisRange(ylo, yhi, "Y");

  if(options.substrate) {
    hSubstrateDens->SetLineColor(kOrange+3); //Brown
    hSubstrateDens->SetLineWidth(2);

    hSubstrateDens->SetAxisRange(xlo, xhi, "X");
    hSubstrateDens->SetAxisRange(ylo, yhi, "Y");
  }
}

void DensFigure::setLegendStyle() {
  hLiquidDens->SetStats(0);
  if(options.substrate) {
    hSubstrateDens->SetStats(0);
  }
}

void DensFigure::setStyle() {
  setLineStyle();
  setLegendStyle();
}

void DensFigure::setValues() {
  if(options.monolayer) {
    monoLimits[0] = dropletPtr->monolayer.zlim[0];
    monoLimits[1] = dropletPtr->monolayer.zlim[1];
  }
}

void DensFigure::drawHists() {
  hLiquidDens->Draw("hist");
  if(options.substrate) {
    hSubstrateDens->Draw("hist same"); //Same canvas
  }
}

void DensFigure::drawLines() {
  if(options.monolayer) {
    monoHiLineDens->SetX1(monoLimits[1]);
    monoHiLineDens->SetX2(monoLimits[1]);
    monoHiLineDens->Draw("same");
    monoLoLineDens->SetX1(monoLimits[0]);
    monoLoLineDens->SetX2(monoLimits[0]);
    monoLoLineDens->Draw("same");
  }
}

void DensFigure::drawLegend() {
  legend->Draw("same");
}

void DensFigure::draw() {
  cout << "setValues" << endl;
  setValues();
  cout << "cd" << endl;
  canvas->cd();
  cout << "setStyle" << endl;
  setStyle();
  cout << "setTitle" << endl;
  setTitle();
  cout << "setAxisLabels" << endl;
  setAxisLabels();
  cout << "drawHists" << endl;
  drawHists();
  cout << "drawLines" << endl;
  drawLines();
  cout << "drawLegend" << endl;
  drawLegend();
  cout << "update canvas" << endl;
  canvas->Update();
  cout << "done drawing" << endl;
}

TanhFigure::TanhFigure(const string _title, TanhFit &tanhFit) : Figure(_title, *tanhFit.simDataPtr) {
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
  legend->AddEntry(tanhFitPtr->hTanh, "data", "l");
  legend->AddEntry(tanhFitPtr->fTanh, "fit", "l");
}

void TanhFigure::setAxisLabels() {
  // col fit
  if(strcmp(rowOrCol, "col") == 0) {
    tanhFitPtr->hTanh->SetXTitle("z (#AA)");
  }
  // row fit
  else {
    tanhFitPtr->hTanh->SetXTitle("r (#AA)");
  }
  tanhFitPtr->hTanh->SetYTitle("mass density (g/cc)");
  tanhFitPtr->hTanh->GetXaxis()->CenterTitle();
  tanhFitPtr->hTanh->GetYaxis()->CenterTitle();
}

void TanhFigure::updateRowCol() {
  strcpy(rowOrCol, tanhFitPtr->rowOrCol);
  rowColNum = tanhFitPtr->rowColNum;
}

void TanhFigure::setTitle() {
  char title[256];
  // Note timestep in title
  sprintf(title, "%s %s %d, t=%09d", simDataPtr->options.dumpfile.data(), rowOrCol, rowColNum, simDataPtr->framePtr->time);
  tanhFitPtr->hTanh->SetTitle(title);
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
  x0Line->SetLineColor(kGreen);
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
  updateRowCol();
  setAxisLabels();
  setTitle();
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
