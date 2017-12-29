#include "Readers.h"

/////////////
// Readers //
/////////////

void InputStream::verifyStream() {
  //Check files
  if (stream.good())
    cout << "Successfully opened " << filename << endl;
  else
    {
      cout << "Failed to open " << filename << endl;
      exit(1);
    }
}

InputStream::~InputStream() {
  stream.close();
}

InputStream::open(char* _filename) {
  filename = _filename;
  stream.open(filename);
  verifyStream();
}

InputStream::InputStream(char* _filename) {
  open(filename);
}


void HeaderReader::setContext(InputStream* _inputStreamPtr, Timestep* _timestepPtr, int* lineNumPtr) {
  inputStreamPtr = _inputStreamPtr;
  timestepPtr = _timestepPtr;
  lineNumPtr = lineNumPtr;
}

void HeaderReader::readHeader() {
  string junk;
  inputStreamPtr->stream >> junk >> timestepPtr->time;
  inputStreamPtr->stream.ignore();
  *lineNumPtr++;

  cout << endl;
  cout << "Timestep " << timestepPtr->time << " @ line " << *lineNumPtr-1 << endl;
}


void LineReader::setContext(InputStream* _inputStreamPtr, int* _atomNumPtr, int* _lineNumPtr) {
  inputStreamPtr = _inputStreamPtr;
  atomNumPtr = _atomNumPtr;
  lineNumPtr = _lineNumPtr;
}

//Split a string into a string vector of words
vector<string> LineReader::strSplit(string str)
{
  int len=str.length();
  stringstream ss(str);
  int numWords=1;
  bool foundSpace=false;

  //Count number of words
  for(int ch=0;ch<len;ch++)
    {
      if(isspace(str[ch]))
        foundSpace=true;
      if(!isspace(str[ch])&&foundSpace)
        {
          numWords++;
          foundSpace=false;
        }
    }

  //Allocate array
  vector<string> arr(numWords);

  //Read string into array
  for(int i=0;i<len;i++)
    ss >> arr[i];

  return arr;
}

//Get index and coordinates from string array of words
void LineReader::strToData(double* coords,double* velocities,double& dipole,string line)
{
  string str;

  int index;
  //Split string into array of words
  vector<string> strArr=strSplit(line);

  //Get index
  str=strArr[0];
  index=atoi(str.data());

  //Save values
  //string -> cstring -> double
  for(int i=0;i<3;i++)
    {
      //Get coordinates from the second, third, and fourth elements in the array
      str=strArr[1+i];
      *coords++=atof(str.data());

      //Get velocities from the sixth, seventh, and eighth elements in the array
      str=strArr[5+i];
      *velocities++=atof(str.data());
    }

  //Save dipole moment cos(theta)
  dipole=atof(strArr[10].data());
}

void LineReader::readLine() {
  double coords[3], velocities[3];
  double dipole;
  string line;

  getline(inputStreamPtr->stream,line);
  *lineNumPtr++;

  //Save data
  strToData(coords,velocities,dipole,line);
  atom.x = coords[0];
  atom.y = coords[1];
  atom.z = coords[2];

  atom.vx = velocities[0];
  atom.vy = velocities[1];
  atom.vz = velocities[2];

  atom.cosTheta = dipole;
  *atomNumPtr++;

  // READ THE LINE AND STORE DATA IN atom
}

void TimestepReader::setContext(InputStream* _inputStreamPtr, AtomArray* atomArrayPtr, Timestep* _timestepPtr, SimData* _simDataPtr) {
  inputStreamPtr = _inputStreamPtr;
  lineReader.setContext(inputStreamPtr, &atomNum, &lineNum);
  headerReader.setContext(inputStreamPtr, timestepPtr, &lineNum);
  timestepPtr = _timestepPtr;
  simDataPTr = _simDataPtr;
}

void TimestepReader::resetAtomCounter() {
  atomNum = 0;
}

void TimestepReader::readTimestep() {
  resetAtomCounter();

  cout << "Step # " << timestepPtr->time << endl;
  headerReader.readHeader();
  for(int i=0; i<atomArrayPtr->numAtoms; i++) {
    lineReader.readLine();
    atomArrayPtr->setAtom(i, lineReader.atom);
  }
}

void FrameReader::openStream(CommandLineParser commandLineParser) {
  inputStream.open(commandLineParser.inLoc);
}

void FrameReader::setContext(CommandLineParser commandLineParser, AtomArray* _atomArrayPtr, Timestep* _timestepPtr, SimData* _simDataPtr) {
  atomArrayPtr = _atomArrayPtr;
  timestepPtr = _timestepPtr;
  simDataPtr = _simDataPtr;
  commandLineParser = _commandLineParser

  stepsPerFrame = simDataPtr->stepsPerFrame;
  openStream(commandLineParser.inLoc);

  timestepReader.setContext(&inputStream, atomArrayPtr, timestepPtr, simDataPtr);
}

FrameReader::FrameReader(CommandLineParser commandLineParser, AtomArray* _atomArrayPtr, Timestep* _timestepPtr, SimData* _simDataPtr) {
  openStream(commandLineParser);
  setContext(_inputStreamPtr, _atomArrayPtr, _timestepPtr);
}

void FrameReader::readFrame() {
  for(int n; n<stepsPerFrame; n++) {
    timestepReader.readTimestep();
  }
}

