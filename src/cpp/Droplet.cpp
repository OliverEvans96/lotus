#include "Droplet.h"

////////////////////////
// Droplet Components //
////////////////////////

Monolayer::Monolayer() {}
Monolayer::~Monolayer() {
  deleteHist();
}

void Monolayer::setContext(Options _options, SimData *_simDataPtr, AtomArray *_atomArrayPtr) {
  options = _options;
  simDataPtr = _simDataPtr;
  atomArrayPtr = _atomArrayPtr;

  tanhFit.setContext(*simDataPtr);
}

void Monolayer::createHist(Grid &grid) {
  hMono = new TH1D("hMono", "hMono", grid.nr, grid.rVals);
  histCreated = true;
  tanhFit.setHist(hMono);
}

void Monolayer::deleteHist() {
  if(histCreated) {
    delete hMono;
  }
}

void Monolayer::calculateRadius() {
  tanhFit.setFitType("mono");
  tanhFit.setFitNum(0);
  tanhFit.solve();
  if(tanhFit.good()) {
    radius = tanhFit.getBoundary();
  }
  else {
    cout << "MONOLAYER RADIUS FIT FAILED" << endl;
  }
}

void Monolayer::reset() {
  hMono->Reset();
}

void Monolayer::fillOne(Atom &atom) {
  double mass;
  mass = simDataPtr->masses[atom.type];
  atom.calculateRadius();
  // TODO: Make this more efficient by calculating
  // only the necessary "radius" (y or r).
  if(options.geometry == "spherical") {
    hMono->Fill(atom.r, mass);
  }
  else if(options.geometry == "cylindrical") {
    hMono->Fill(atom.y, mass);
  }
}

void Monolayer::fill(AtomArray &atoms) {
  Atom atom;
  reset();
  // Add to monolayer
  for(int stepInFrame=0; stepInFrame<simDataPtr->framePtr->stepsThisFrame; stepInFrame++) {
    for(int atomNum=0; atomNum<simDataPtr->numAtoms; atomNum++) {
      atoms.getAtom(atomNum, stepInFrame, atom);
      if(isIn(atom.type, simDataPtr->liquidTypes)) {
        if(inMonolayer(atom)) {
          fillOne(atom);
        }
      }
    }
  }

  convertUnits();
}

void Monolayer::convertUnits() {
  // Must be called after monoLimits are found.
  double dv;

  // Divide by number of steps per frame
  hMono->Scale(1.0/simDataPtr->framePtr->stepsThisFrame);
  // Calculate monolayer bin volume
  // options.dv is for hDroplet, which uses options.dz.
  dv = options.dv * (zlim[1] - zlim[0])/options.dz;
  // Divide by volume to get density
  hMono->Scale(1.0/dv);
  // Convert units from amu/AA to g/cc
  hMono->Scale(NANO_DENS_TO_MACRO);
  // Reset error created by scaling
  for(int i=1; i<=hMono->GetNbinsX(); i++) {
    hMono->SetBinError(i, 0.0);
  }
}

// Whether an atom is in the monolayer
bool Monolayer::inMonolayer(Atom &atom) {
  return ((zlim[0] < atom.z) && (atom.z < zlim[1]));
}

//Keep track of which atoms join the monolayer, and save their radius scaled by the base radius to the vector rScaledJoin
int Monolayer::monoFlux(vector<double> r,vector<double> z,double* monoLimits,double baseRadius,TH1D* rScaledJoin,TH1D* rScaledLeave,int &nMono)
{
    //Number of atoms in the simulation
    int numAtoms=z.size();

    //List of atoms currently in the monolayer
    vector<int> monoList;
    //List of atoms in the monolayer last timestep
    static    vector<int> prevMonoList;

    //Calculate flux since last step
    int flux = 0;

    //Count how many molecules should be in the monolayer based on the flux.
    static int expected = 0;

    //Save ids of atoms below the bulk-mono threshold
    for(int i=0;i<numAtoms;i++)
    {
        //if(i%(numAtoms/20)==0)
        //    cout << "#" << i << ": " << z[i] << endl;
        if((monoLimits[0]<=z[i]) && (z[i]<=monoLimits[1]))
            monoList.push_back(i);
    }

    //Number of atoms in the monolayer
    int numMonoAtoms=monoList.size();
    //Number of atoms in the monolayer last step
    int numPrevMonoAtoms=prevMonoList.size();

    //Save number of molecules in monolayer
    nMono = numMonoAtoms;


    cout << endl;
    cout << "Prev monolayer: " << numPrevMonoAtoms << endl;
    cout << "This monolayer: " << numMonoAtoms << endl;
    cout << endl;

    //For each atom in the monolayer, check whether it was there last step
    for(int i=0;i<numMonoAtoms;i++)
    {
        //Join
        //If it's new to the monolayer, save it's scaled radial position
        if( !isIn(monoList[i],prevMonoList) )
        {
            flux++;
            rScaledJoin->Fill(r[i]/baseRadius);
        }
    }

    //For each atom in the monolayer last time, check whether it's still there
    for(int i=0;i<numPrevMonoAtoms;i++)
    {
        //Leave
        //If has just left the monolayer, save it's scaled radial position
        if( !isIn(prevMonoList[i],monoList) )
        {
            flux--;
            rScaledLeave->Fill(r[i]/baseRadius);
        }
    }

    //Check whether the expected value matches the counted value
    /*
    expected += flux;
    cout << "Expected: " << expected;
    cout << " Counted: " << numMonoAtoms;
    if (expected == numMonoAtoms)
        cout << " OK!" << endl;
    else
        cout << " NOOOOOOOOOOOO!" << endl;
    */

    //Save monoList for next time
    prevMonoList=monoList;

    return flux;
}

//Find z limits on monolayer
void Monolayer::findMonoLimits(TH1D *hLiquidDens,double *monoLimits)
{
    int n = hLiquidDens->GetNbinsX();
    bool foundPeak=false;
    bool foundDip=false;
    double dens=0;

    for(int i=1;i<=n;i++)
    {
        dens=hLiquidDens->GetBinContent(i);

        //Find first peak above 1 - this is the beginning of the monolayer
        if(!foundPeak&&dens>1)
        {
            foundPeak=true;
            cout << "Beginning monolayer at " << hLiquidDens->GetBinLowEdge(i) << endl;
            monoLimits[0]=hLiquidDens->GetBinLowEdge(i);
        }

        //if(foundPeak)
        //    cout << "dens=" << dens << " at " << hLiquidDens->GetBinLowEdge(i) << endl;

        if(foundPeak&&dens>1)
            cout << "Still in monolayer at " << hLiquidDens->GetBinLowEdge(i) << endl;

        //Find first drop below 1 after peak - this is the end of the monolayer
        if(foundPeak&&dens<=1)
        {
            foundDip=true;
            //cout << "Monolayer ends at " << hLiquidDens->GetBinLowEdge(i) << endl;
            monoLimits[1]=hLiquidDens->GetBinLowEdge(i);
            break;
        }
    }

    if(foundDip)
        cout << "Found monolayer limits: " << monoLimits[0] << " " << monoLimits[1] << endl;
    else
        cout << "FAILED to locate monolayer." << endl;

}

CircularBulk::CircularBulk() {
  gCirclePoints = new TGraph();
  numPoints = 0;
  headers[0] = "x";
  headers[1] = "y";
}

CircularBulk::~CircularBulk() {
  delete gCirclePoints;
}

void CircularBulk::setContext(Options _options, SimData *_simDataPtr, AtomArray *_atomArrayPtr) {
  options = _options;
  simDataPtr = _simDataPtr;
  atomArrayPtr = _atomArrayPtr;
  tanhFit.setContext(*simDataPtr);
  circle.setContext(*simDataPtr, gCirclePoints);
}

void CircularBulk::setHist(TH2D *_hDroplet) {
  hDroplet = _hDroplet;
}

bool CircularBulk::pointOk(double r, double z) {
  bool ok = true;

  if(z < simDataPtr->monoTop) {
    ok = false;
  }

  return ok;
}

void CircularBulk::saveBoundaryPoints() {
  numPoints = gCirclePoints->GetN();
  boundaryPointsArray[0] = gCirclePoints->GetX();
  boundaryPointsArray[1] = gCirclePoints->GetY();
}

void CircularBulk::findBoundaryPoints() {
  int nx, ny;
  double x, y;
  int pointNum;

  nx = hDroplet->GetNbinsX();
  ny = hDroplet->GetNbinsY();

  gCirclePoints->Set(0);
  pointNum = 0;

  // Row fits (start at first row above monolayer)
  tanhFit.setFitType("row");
  for(int j=firstBulkBin; j<=ny; j++) {
    // cout << "j = " << j << endl;
    tanhFit.setHist(hDroplet->ProjectionX("px", j, j));
    tanhFit.setFitNum(j);
    // cout << "solve" << endl;
    tanhFit.solve();
    // cout << "verify" << endl;
    // Solve fails if droplet is not present in this row.
    if(tanhFit.good()) {
      x = tanhFit.getBoundary();
      y = hDroplet->GetYaxis()->GetBinCenter(j);
      if(pointOk(x, y)) {
        gCirclePoints->SetPoint(pointNum++, x, y);
      }
    }
  }

  // Column fits
  tanhFit.setFitType("col");
  for(int i=1; i<=nx; i++) {
    tanhFit.setHist(hDroplet->ProjectionY("py", i, i));
    tanhFit.setFitNum(i);
    // Solve fails if droplet is not present in this row.
    tanhFit.solve();
    if(tanhFit.good()) {
      x = hDroplet->GetXaxis()->GetBinCenter(i);
      y = tanhFit.getBoundary();
      if(pointOk(x, y)) {
        gCirclePoints->SetPoint(pointNum++, x, y);
      }
    }
  }

  saveBoundaryPoints();
}

//Given a graph of points, a maximum x value, and the y coordinate of the interface with the substrate, fit a circle to the points and find the intersection of the circle with the substrate interface. The result is the bulk-monolayer interface
double CircularBulk::fitCircle() {
  gCirclePoints->Sort();

  circle.fit();
  circle.Print();

  // Update points stored for writing
  // since fitting may discard some
  saveBoundaryPoints();

  //Get chi2 scaled by number of points in fit
  return circle.GetChi2s();
}

Droplet::Droplet(AtomArray &atomArray) {
  // TODO: Determine radius from options?
  rDensCyl = 10;

  setContext(atomArray);
  createHists();
}

Droplet::~Droplet() {
  delete hDroplet;
}

void Droplet::setContext(AtomArray &atomArray) {
  atomArrayPtr = &atomArray;
  simDataPtr = atomArrayPtr->simDataPtr;
  options = simDataPtr->options;

  bulk.setContext(options, simDataPtr, atomArrayPtr);
  monolayer.setContext(options, simDataPtr, atomArrayPtr);
}

void Droplet::fillOne(Atom &atom) {
  double mass, z;

  mass = simDataPtr->masses[atom.type];
  // z relative to top of substrate
  z = atom.z - simDataPtr->substrateTop;
  // cout << "atom.z = " << atom.z << endl;
  // cout << "substrateTop = " << simDataPtr->substrateTop << endl;
  // cout << "z = " << z << endl;

  // 2D Density
  hDroplet->Fill(atom.r, z, mass);
  // 1D Density
  if(atom.r < rDensCyl) {
    hLiquidDens->Fill(atom.z, mass);
  }

  // cout << "Filling type " << atom.type << " (droplet): " << mass << " @ " << atom.z << endl;
}

void Droplet::fill(AtomArray &atoms) {
  Atom atom;
  reset();
  for(int stepInFrame=0; stepInFrame<simDataPtr->framePtr->stepsThisFrame; stepInFrame++) {
    for(int atomNum=0; atomNum<simDataPtr->numAtoms; atomNum++) {
      atoms.getAtom(atomNum, stepInFrame, atom);
      // If liquid
      if(isIn(atom.type, simDataPtr->liquidTypes)) {
        atom.calculateRadius();
        fillOne(atom);
      }
    }
  }

  convertUnits();
}

void Droplet::findMonolayer() {
  monolayer.findMonoLimits(hLiquidDens, monolayer.zlim);
  // Determine first bulk row to use for row fits
  simDataPtr->monoTop = monolayer.zlim[1]-simDataPtr->substrateTop;
  bulk.firstBulkBin = hDroplet->GetYaxis()->FindBin(simDataPtr->monoTop)+1;
  cout << "monolayer.zlim: " << monolayer.zlim[0] << " " << monolayer.zlim[1] << endl;
}

void Droplet::reset() {
  hDroplet->Reset();
  hLiquidDens->Reset();
}

void Droplet::convertUnits() {
  // Divide by number of steps per frame
  hDroplet->Scale(1.0/simDataPtr->framePtr->stepsThisFrame);
  hLiquidDens->Scale(1.0/simDataPtr->framePtr->stepsThisFrame);
  // Divide by volume to get density
  hDroplet->Scale(1.0/options.dv);
  hLiquidDens->Scale(1.0/(PI*options.dz*square(rDensCyl)));
  // Convert units from amu/AA^3 to g/cc
  hDroplet->Scale(NANO_DENS_TO_MACRO);
  hLiquidDens->Scale(NANO_DENS_TO_MACRO);
  // Reset error created by scaling
  for(int i=1; i<=hLiquidDens->GetNbinsX(); i++) {
    hLiquidDens->SetBinError(i, 0.0);
  }
  for(int i=1; i<=hDroplet->GetNbinsY(); i++) {
    for(int j=1; j<=hDroplet->GetNbinsX(); j++) {
      hDroplet->SetBinError(i, j, 0.0);
    }
  }
}

void Droplet::createHists() {
  double zlo, zhi;
  double rlo, rhi;
  // Set limits on histograms
  Grid grid;
  int nz;
  double dz = options.dz;
  double dv = options.dv;
  zlo = 0;
  zhi = options.plot_zmax;
  rlo = 0;
  rhi = options.plot_rmax;

  grid.setBounds(zlo, zhi, rhi);
  grid.setSpacing(dz, dv);
  grid.init();

  if(options.verbose) {
    cout << "Grid:" << endl;
    cout << "dz = " << dz << endl;
    cout << "dv = " << dv << endl;
    cout << "nr = " << grid.nr << endl;
    cout << "nz = " << grid.nz << endl;
    cout << "zlo = " << grid.zlo << endl;
    cout << "zhi = " << grid.zhi << endl;
    cout << "rVals @ " << grid.rVals << " = ";
    for(int i=0; i<grid.nr; i++) {
      cout << grid.rVals[i] << ", ";
    }
    cout << endl;
  }

  hDroplet = new TH2D("hDroplet","hDroplet",grid.nr,grid.rVals,grid.nz,grid.zlo,grid.zhi);
  hDroplet->SetStats(0);
  hDroplet->SetMinimum(0);
  hDroplet->SetMaximum(2.0);
  bulk.setHist(hDroplet);

  monolayer.createHist(grid);

  // TODO: This is pretty messy.
  // Want to use actual coords for this one
  // Rather than shifted z coords as above
  zlo = simDataPtr->simBounds.zlo;
  zhi = simDataPtr->simBounds.zhi;
  nz = (int) ceil((zhi - zlo)/dz);
  // If dz doesn't evenly divide zhi-zlo, shift zhi up slightly.
  zhi = zlo + nz*dz;
  hLiquidDens = new TH1D("hLiquidDens", "hLiquidDens", nz, zlo, zhi);
}

double Droplet::getMass() {
  // Total mass in amu
  double mass;

  // Integrate over whole cylindrical volume
  // All bins have equal volume by design
  // This quantity has units g/cm^3*AA^3
  mass = hDroplet->Integral() * options.dv;
  // Convert units to amu
  mass /= NANO_DENS_TO_MACRO;
  return mass;
}

double Droplet::getMass1D() {
  double mass;
  // TODO: Check scaling on hLiquidDens.
  // TODO: This is actually just mass in cylinder, not total
  mass = hLiquidDens->Integral() * options.dz * PI * square(rDensCyl);
  mass /= NANO_DENS_TO_MACRO;
  return mass;
}

void Droplet::dropletCalculations() {
  // TODO: Set boundary points from options
  double rBulkMax = 200.0;
  double chi2;

  // Get circle
  bulk.findBoundaryPoints();
  if(options.fitCircle) {
    chi2 = bulk.fitCircle();
    bulk.radius = bulk.circle.Intersect(simDataPtr->monoTop);
    bulk.contactAngle = bulk.circle.GetContactAngle();
    bulk.height = bulk.circle.GetHeight();
  }

  cout << "monoTop = " << simDataPtr->monoTop << endl;
  cout << "radius = " << bulk.radius << endl;
  cout << "contactAngle = " << bulk.contactAngle << endl;
  cout << "height = " << bulk.height << endl;
}
