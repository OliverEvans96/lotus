#include "MDBase.h"

//Count the number of water atoms in the first timestep
int SimData::countAtoms(ifstream &inFile)
{
  //Count number of atoms
  bool countFlag=true;
  string line;
  int numAtoms=0;

  //Ignore the first 3 lines
  for(int i=0;i<3;i++) inFile.ignore(256,'\n');

  while(countFlag)
    {
      getline(inFile,line);

      //Count until reaching a line containing "TIMESTEP"
      if(line.find("TIMESTEP")!=string::npos||inFile.eof())
        {
          countFlag=false;
          numAtoms-=1; //Account for the blank line
        }
      else
        numAtoms++;
    }

  //Unset eof flag (if set) & return to beginning of file
  inFile.clear();
  inFile.seekg(0,ios::beg);

  cout << "Counted " << numAtoms << " atoms." << endl;

  return numAtoms;
}

//Count the number of timesteps
int SimData::countSteps(ifstream &inFile)
{
  string line;
  int numSteps=0;
  int lineNum=0;

  //Count number of timesteps
  while(getline(inFile,line))
    {
      if(line.find("TIMESTEP")!=string::npos)
        numSteps++;
    }

  //Unset eof flag (if set) & return to beginning of file
  inFile.clear();
  inFile.seekg(0,ios::beg);

  cout << "Counted " << numSteps << " timesteps." << endl;
  return numSteps;
}

SimData::SimData(Options options) {
  inputStream.open(options.inLoc);
  numAtoms = countAtoms(inputStream.stream);
  numSteps = countSteps(inputStream.stream);
}

void SimData::setStepsPerFrame(int _stepsPerFrame) {
  stepsPerFrame = _stepsPerFrame;
  lastFrame.setSteps(numSteps, stepsPerFrame);

  //Number of frames (collections of timesteps)
  //If not divisible, then the last frame will have more than the rest
  numFrames = (int) floor(numSteps/stepsPerFrame);
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
  zhi = _zhi;
  rhi_preround = _rhi;
}

void Grid::setSpacing(double _dz, double _dv) {
  dz = _dz;
  dv = _dv;
}

void Grid::calculateVolumeLimits() {
  // Volume of full cylinder containing outer bin.
  vhi_preround = dz * PI * rhi_preround * rhi_preround;
  nv = (int) round(vhi_preround / dv);
  nr = nv;
  vhi = nv * dv;
  rhi = sqrt(vhi / (PI*dz));

  nz = (int) round((zhi_preround-zlo) / dz);
  zhi = zlo + nz * dz;
}

void Grid::allocateBins() {
  rVals = new double[nv];
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

