#include "MDBase.h"

SimData::SimData(Options options) {
  setOptions(options);

  // Default values
  substrateTop = 0;
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

void SimData::setOptions(Options _options) {
  options = _options;
  liquidTypes = options.liquidTypes;
  solidTypes = options.solidTypes;
  stepsPerFrame = options.stepsPerFrame;

  // Optionally specify number of atoms by hand
  // e.g. if not all atoms are dumped, causing
  // dumpfile and datafile to have different #s.
  if(options.numAtoms > 0) {
    numAtoms = options.numAtoms;
  }
}

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

void SimData::setNumSteps(int _numSteps) {
  numSteps = _numSteps;
  setStepsPerFrame(options.stepsPerFrame);
}

// void Grid::mmmmm() {
//   // Replace with Grid object ////////////////////////////////
//   double dz=1.5;                                            //
//   double zlo=10;                                            //
//   double zhi=142; //142 for 60A, 121 for others             //
//   int nz=(int) round((zhi-zlo)/dz);                         //
//                                                             //
//   //Check that dz, zlo, and zhi are compatible              //
//   if( !(fmod((zhi-zlo),dz)==0) )                            //
//   {                                                         //
//     cout << "(zhi-zlo) is not divisible by dz" << endl;     //
//     return 1;                                               //
//   }                                                         //
//                                                             //
//   cout << "nz=" << nz << endl;                              //
//                                                             //
//   const double dA=500;                                      //
//   const double alo=0;                                       //
//   const double ahi=125500;                                  //
//   int nA=(int) round((ahi-alo)/dA);                         //
//                                                             //
//   double dr=dz;                                             //
//   double rlo=sqrt(alo/PI);                                  //
//   double rhi=sqrt(ahi/PI);                                  //
//   int nr=(int) round((rhi-rlo)/dr);                         //
//                                                             //
//   double vlo=0;                                             //
//   double vhi=1e7;                                           //
//   int nV=(int) round((vhi-vlo)/dA);                         //
//                                                             //
//   double dp=dr;                                             //
//   double plo=pow((3*vlo/(4*PI)),(1.0/3));                   //
//   double phi=pow((3*vhi/(4*PI)),(1.0/3));                   //
//   int  np=(int) round((phi-plo)/dp);                        //
//                                                             //
//   //Bin volume                                              //
//   double dV=dA*dz;                                          //
//                                                             //
//   //Check that dA, alo, and ahi are compatible              //
//   if( !(fmod((ahi-alo),dA)==0) )                            //
//   {                                                         //
//     cout << "(Ahi-Alo) is not divisible by dA" << endl;     //
//     return 1;                                               //
//   }                                                         //
//                                                             //
//   cout << "nA=" << nA << endl;                              //
//                                                             //
//   ////////////////////////////////////////////////////////////
//
//
//   //Values to use for equal-area bins
//   double aVals[nA];
//   double vVals[nA];
//   double rVals[nA];
//   double pVals[nA];
//
// }

void Grid::setBounds(double _zlo, double _zhi, double _rhi) {
  zlo = _zlo;
  zhi_preround = _zhi;
  rhi_preround = _rhi;
}

void Grid::setSpacing(double _dz, double _dv) {
  dz = _dz;
  dv = _dv;
}

void Grid::calculateVolumeLimits() {
  // Volume of full cylinder containing outer bin.
  // Increase upper limits if total width is not divisible
  // by required spacing (e.g. (zhi-zlo)%dz != 0)
  vhi_preround = dz * PI * rhi_preround * rhi_preround;
  nv = (int) ceil(vhi_preround / dv);
  nr = nv;
  vhi = nv * dv;
  rhi = sqrt(vhi / (PI*dz));
  nz = (int) ceil((zhi_preround-zlo) / dz);
  zhi = zlo + nz * dz;
}

void Grid::allocateBins() {
  // Includes both endpoints
  // https://root.cern.ch/root/roottalk/roottalk10/1170.html
  rVals = new double[nv+1];
}

void Grid::calculateBins() {
  for(int i=0;i<=nv;i++) {
    rVals[i]=sqrt(i*dv/(PI*dz));
    // pVals[i]=pow(((3*vVals[i])/(4*PI)),1.0/3);
  }
}

void Grid::init() {
  calculateVolumeLimits();
  allocateBins();
  calculateBins();
}

void Grid::deallocateBins() {
  delete [] rVals;
}

Grid::~Grid() {
  deallocateBins();
}
