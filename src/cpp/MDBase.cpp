#include "MDBase.h"

SimData::SimData(Options options) {
  setOptions(options);
}

SimData::~SimData() {
  deleteWaterBonds();
}

void SimData::deleteWaterBonds() {
  map<int, int*>::iterator it;
  for(it=waterBonds.begin(); it!=waterBonds.end(); it++) {
    delete it->second;
  }
}

/**
   Copy a few important options directly to this object.

   @see #liquidTypes
   @see #solidTypes
   @see #stepsPerFrame
   @see #monoTop
   @see #substrateTop
   @see #numAtoms
*/
void SimData::setOptions(Options _options) {
  options = _options;
  liquidTypes = options.liquidTypes;
  solidTypes = options.solidTypes;
  stepsPerFrame = options.stepsPerFrame;
  monoTop = options.monoTop;
  substrateTop = options.substrateTop;

  // Optionally specify number of atoms by hand
  // e.g. if not all atoms are dumped, causing
  // dumpfile and datafile to have different #s.
  if(options.numAtoms > 0) {
    numAtoms = options.numAtoms;
  }
}

/**
   Set the number of steps per frame and determine
   how many steps are in the last frame.
   Must be called after #numSteps is defined.
   This function is called by #setNumSteps.
   @see LastFrame::setSteps
*/
void SimData::setStepsPerFrame(int _stepsPerFrame) {
  stepsPerFrame = _stepsPerFrame;
  lastFrame.setSteps(numSteps, stepsPerFrame);
  if(options.verbose) {
    cout << "Last frame:" << endl;
    cout << "extraSteps = " << lastFrame.extraSteps << endl;
    cout << "numSteps = " << lastFrame.numSteps << endl;
    cout << "frameNum = " << lastFrame.frameNum << endl;
  }

  //Number of frames (collections of timesteps)
  //If not divisible, then the last frame will have more than the rest
  numFrames = (int) floor(numSteps/stepsPerFrame);
}

/**
   Set numsteps and call #setStepsPerFrame
   using value from Options::stepsPerFrame.
*/
void SimData::setNumSteps(int _numSteps) {
  numSteps = _numSteps;
  setStepsPerFrame(options.stepsPerFrame);
}

void Grid::setBounds(double _zlo, double _zhi, double _rhi) {
  zlo = _zlo;
  zhi_preround = _zhi;
  rhi_preround = _rhi;
}

void Grid::setSpacing(double _dz, double _dv) {
  dz = _dz;
  dv = _dv;
}

/**
   Increase upper limits if total width is not divisible
   by required spacing (i.e. (zhi-zlo)%dz != 0)
*/
void Grid::calculateVolumeLimits() {
  // Volume of full cylinder containing outer bin.
  vhi_preround = dz * PI * rhi_preround * rhi_preround;
  nv = (int) ceil(vhi_preround / dv);
  nr = nv;
  vhi = nv * dv;
  rhi = sqrt(vhi / (PI*dz));
  nz = (int) ceil((zhi_preround-zlo) / dz);
  zhi = zlo + nz * dz;
}

/**
   Includes both endpoints
   https://root.cern.ch/root/roottalk/roottalk10/1170.html
*/
void Grid::allocateBins() {
  rVals = new double[nv+1];
  allocated = true;
}

/**
   @f$r_i = \sqrt{\frac{i\,dv}{\pi\,dz}} @f$
*/
void Grid::calculateBins() {
  for(int i=0;i<=nv;i++) {
    rVals[i]=sqrt(i*dv/(PI*dz));
  }
}

/**
   Calls #calculateVolumeLimits, #allocateBins, #calculateBins
*/
void Grid::init() {
  calculateVolumeLimits();
  allocateBins();
  calculateBins();
}

void Grid::deallocateBins() {
  if(allocated) {
    delete [] rVals;
  }
}

/**
   Call #deallocateBins
*/
Grid::~Grid() {
  deallocateBins();
}
