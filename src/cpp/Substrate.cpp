#include "Substrate.h"

Substrate::Substrate(AtomArray &atomArray, double dz) {
  setContext(atomArray);
  createHist(dz);
}

Substrate::~Substrate() {
  delete hSubstrateDens;
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
  hSubstrateDens->Fill(atom.z, mass);
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
  hSubstrateDens->Reset();
}

void Substrate::convertUnits() {
  double dx, dy, dz;

  // Divide by number of steps per frame
  hSubstrateDens->Scale(1.0/simDataPtr->stepsPerFrame);
  // Divide by volume to get density
  dx = simDataPtr->simBounds.xhi - simDataPtr->simBounds.xlo;
  dy = simDataPtr->simBounds.yhi - simDataPtr->simBounds.ylo;
  dz = hSubstrateDens->GetXaxis()->GetBinWidth(0);
  hSubstrateDens->Scale(1.0/(dx*dy*dz));
  // Convert units from amu/AA^3 to g/cc
  hSubstrateDens->Scale(NANO_DENS_TO_MACRO);
}

void Substrate::createHist(double dz) {
  double zlo, zhi;
  int nz;
  zlo = simDataPtr->simBounds.zlo;
  zhi = simDataPtr->simBounds.zhi;
  nz = (int) ceil((zhi - zlo)/dz);
  // If dz doesn't evenly divide zhi-zlo, shift zhi up slightly.
  zhi = zlo + nz*dz;

  hSubstrateDens = new TH1D("hSubstrateDens", "hSubstrateDens", nz, zlo, zhi);
}

void Substrate::findLimits() {
  bool foundBottom = false;
  // Cutoff density (g/cc)
  double cutoff = 0.05;
  double dens;

  // Require this many bins under the cutoff
  // in order to consider the top of the substrate found.
  int numRequired = 3;
  int numFound = 0;
  // ROOT hists are 1-indexed
  for(int i=1; i<=hSubstrateDens->GetNbinsX(); i++) {
    dens = hSubstrateDens->GetBinContent(i);
    if(!foundBottom && dens>cutoff) {
      zlim[0] = hSubstrateDens->GetBinLowEdge(i);
      printf("Found bottom @ bin %d (z=%.2f)\n", i, zlim[0]);
      foundBottom = true;
    }
    if((foundBottom) && dens<cutoff) {
      numFound++;
    }
    if(numFound == numRequired) {
      // If we have seen enough bins below cutoff,
      // then call the first one the top.
      zlim[1] = hSubstrateDens->GetBinLowEdge(i-(numRequired-1));
      printf("Found top @ bin %d (z=%.2f)\n", i, zlim[1]);
      break;
    }
  }
  if(options.verbose) {
    printf("Substrate limits: (%.2f, %.2f)\n", zlim[0], zlim[1]);
  }

  // Save to simData
  cout << "Setting substrateTop = " << zlim[1] << endl;
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
  mass = hSubstrateDens->Integral("width") * dx * dy;
  // Convert units to amu
  mass /= NANO_DENS_TO_MACRO;
  return mass;
}

