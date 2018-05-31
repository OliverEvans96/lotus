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
    // If not liquid
    if(!isIn(atoms.type[i], simDataPtr->liquidTypes)) {
      atoms.getAtom(i, atom);
      fillOne(atom);
    }
  }
  convertUnits();
}

void Substrate::reset() {
  hSubs->Reset();
}

void Substrate::convertUnits() {
  double dx, dy, dz;

  // Divide by number of steps per frame
  hSubs->Scale(simDataPtr->stepsPerFrame);
  // Divide by volume to get density
  // TODO (get volume)
  dx = simDataPtr->simBounds.xhi - simDataPtr->simBounds.xlo;
  dy = simDataPtr->simBounds.yhi - simDataPtr->simBounds.ylo;
  dz = hSubs->GetXaxis()->GetBinWidth(0);
  hSubs->Scale(1.0/(dx*dy*dz));
  // Convert units from amu/AA to g/cc
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
  cout << "zlo = " << zlo << endl;
  cout << "zhi = " << zhi << endl;
  cout << "dz = " << dz << endl;
  cout << "nz = " << nz << endl;
  hSubs = new TH1D("Substrate", "Substrate", nz, zlo, zhi);

  hSubs->SetLineColor(kOrange+3); //Brown
  hSubs->SetLineWidth(2);
}

void Substrate::createCanvas(int width, int height) {
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

