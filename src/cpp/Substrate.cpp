#include "Substrate.h"

Substrate::Substrate(AtomArray &atomArray, double dz) {
  setContext(atomArray);
  createHist(dz);
  createCanvas();
}

Substrate::~Substrate() {
  delete hSubs;
  delete cSubs;
}

void Substrate::setContext(AtomArray &atomArray) {
  atomArrayPtr = &atomArray;
  simDataPtr = atomArrayPtr->simDataPtr;
  options = simDataPtr->options;
}

void Substrate::fillOne(Atom &atom) {
  int mass;
  mass = simDataPtr->masses[atom.type];
  hSubs->Fill(atom.z, mass);
  // cout << "Filling type " << atom.type << ": " << mass << " @ " << atom.z << endl;
}

void Substrate::fill(AtomArray &atoms) {
  Atom atom;
  reset();
  for(int i=0; i<simDataPtr->numAtoms; i++) {
    // If solid
    if(isIn(atoms.type[i], simDataPtr->solidTypes)) {
      atoms.getAtom(i, atom);
      fillOne(atom);
    }
  }
  convertUnits();
  findLimits();
}

void Substrate::reset() {
  hSubs->Reset();
}

void Substrate::convertUnits() {
  double dx, dy, dz;

  // Divide by number of steps per frame
  hSubs->Scale(1.0/simDataPtr->stepsPerFrame);
  // Divide by volume to get density
  dx = simDataPtr->simBounds.xhi - simDataPtr->simBounds.xlo;
  dy = simDataPtr->simBounds.yhi - simDataPtr->simBounds.ylo;
  dz = hSubs->GetXaxis()->GetBinWidth(0);
  hSubs->Scale(1.0/(dx*dy*dz));
  // Convert units from amu/AA^3 to g/cc
  hSubs->Scale(NANO_DENS_TO_MACRO);
}

void Substrate::createHist(double dz) {
  double zlo, zhi;
  int nz;
  zlo = simDataPtr->simBounds.zlo;
  zhi = simDataPtr->simBounds.zhi;
  nz = (int) ceil((zhi - zlo)/dz);
  // If dz doesn't evenly divide zhi-zlo, shift zhi up slightly.
  zhi = zlo + nz*dz;

  hSubs = new TH1D("Substrate", "Substrate", nz, zlo, zhi);
  hSubs->SetLineColor(kOrange+3); //Brown
  hSubs->SetLineWidth(2);
  hSubs->SetStats(0);
}

// TODO: Move to Visualize.cpp
void Substrate::createCanvas() {
  int width, height;
  width = options.plot_width;
  height = (int) round(width / options.plot_aspect);
  cSubs = new TCanvas("Substrate", "Substrate", width, height);
}

// TODO: Move to Visualize.cpp
void Substrate::plotDensity(char* filename) {
  stringstream ss;
  ss << options.outLoc << "/" << filename;
  gStyle->SetCanvasPreferGL(true);
  cSubs->cd();
  hSubs->Draw("L");
  cSubs->SaveAs(ss.str().data());
  if(options.verbose) {
    cout << "Saving density hist @ '" << ss.str().data() << "'" << endl;
  }
}

void Substrate::findLimits() {
  bool foundBottom = false;
  bool foundTop = false;
  // Cutoff density (g/cc)
  double cutoff = 0.05;
  double dens;
  // ROOT hists are 1-indexed
  for(int i=1; i<=hSubs->GetNbinsX(); i++) {
    dens = hSubs->GetBinContent(i);
    if(!foundBottom && dens>cutoff) {
      zlim[0] = hSubs->GetBinLowEdge(i);
      printf("Found bottom @ bin %d (z=%.2f)\n", i, zlim[0]);
      foundBottom = true;
    }
    if((foundBottom && !foundTop) && dens<cutoff) {
      zlim[1] = hSubs->GetBinLowEdge(i);
      printf("Found top @ bin %d (z=%.2f)\n", i, zlim[1]);
      foundTop = true;
      break;
    }
  }
  if(options.verbose) {
    printf("Substrate limits: (%.2f, %.2f)\n", zlim[0], zlim[1]);
  }

  // Save to simData
  simDataPtr->substrateTop = zlim[1];
}

double Substrate::getMass() {
  // Total mass in amu
  double dx, dy, mass;
  dx = simDataPtr->simBounds.xhi - simDataPtr->simBounds.xlo;
  dy = simDataPtr->simBounds.yhi - simDataPtr->simBounds.ylo;

  // integral gives area density for whole frame
  // multiply by dx, dy to get mass in g/cm^3*AA^3
  // "width" option does integral instead of sum
  mass = hSubs->Integral("width") * dx * dy;
  // Convert units to amu
  mass /= NANO_DENS_TO_MACRO;
  return mass;
}

