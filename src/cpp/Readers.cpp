#include "Readers.h"

//////////////////
// Input Stream //
//////////////////

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

InputStream::InputStream() {};

InputStream::~InputStream() {
  stream.close();
}

void InputStream::open(string _filename) {
  lineNum = 0;
  filename = _filename;
  stream.open(filename);
  verifyStream();
}

InputStream::InputStream(string _filename) {
  open(_filename);
}

void InputStream::skipLines(int numLines) {
  // Ignore lines from the file (at most 256 characters)
  for(int i=0; i<numLines; i++) {
    stream.ignore(256, '\n');
  }

  // Increment line number counter accordingly
  lineNum += numLines;
}

void InputStream::skipWhitespace() {
  int num;
  num = stream.tellg();
  while(isspace(stream.peek())) {
    stream.get();
  }
  num = stream.tellg();
}

bool InputStream::search(string term) {
  // Move the ifstream cursor to the beginning of
  // the first line containing `term`.
  // (the matching line will not be consumed)
  // Returns whether the line was found

  string line;
  bool found;
  int pos;

  found = false;

  while(stream.good()) {
    // Save beginning of line
    pos = stream.tellg();
    getline(stream, line);
    if(isIn(term, line)) {
      found = true;
      // Rewind to line beginning
      stream.seekg(pos);
      break;
    }
  }

  return found;
}

string InputStream::search(vector<string> terms) {
  // Move the ifstream cursor one line past
  // the first line containing one of the
  // strings in `terms`.
  // (the matching line will be consumed)
  // Returns the line found

  string line;
  string term;
  bool found;
  vector<string>::iterator it;

  term = "";
  found = false;

  while(stream.good()) {
    getline(stream, line);
    for(it=terms.begin(); it<terms.end(); it++) {
      term = *it;
      if(isIn(line, term)) {
        found = true;
        break;
      }
    }
  }

  if(!found) {
    term = "";
  }

  return term;
}

bool InputStream::searchLine(string term) {
  // Move the ifstream cursor one line past
  // the first line exactly equal to `term`.
  // (the matching line will be consumed)
  // Returns whether the line was found

  string line;
  bool found;

  found = false;

  while(stream.good()) {
    getline(stream, line);
    if(line == term) {
      found = true;
      break;
    }
  }

  return found;
}

string InputStream::searchLine(vector<string> terms) {
  // Move the ifstream cursor one line past
  // the first line exactly equal to one of the
  // strings in `terms`.
  // (the matching line will be consumed)
  // Returns the line found

  string line;
  string term;
  bool found;
  vector<string>::iterator it;
  int i = 0;

  term = "";
  found = false;

  while(!found && stream.good()) {
    getline(stream, line);
    i++;
    for(it=terms.begin(); it<terms.end(); it++) {
      term = *it;
      if(line == term) {
        found = true;
        break;
      }
    }
  }

  if(!found) {
    term = "";
  }

  return term;
}

bool InputStream::nextLineBlank() {
  bool blank = false;
  if(stream.peek() == '\n') {
    stream.get();
    blank = (stream.peek() == '\n');
  }
  return blank;
}

string InputStream::peekLine() {
  string line;
  int pos = stream.tellg();
  getline(stream, line);
  stream.seekg(pos);
  return line;
}

//Split a string into a string vector of words
vector<string> strSplit(string str)
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

// TODO: Delete
void strToData(double *coords,double *velocities,double &dipole,string line)
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

/////////////
// Readers //
/////////////

void HeaderReader::setContext(InputStream* _inputStreamPtr, Timestep* _timestepPtr, int* _lineNumPtr) {
  inputStreamPtr = _inputStreamPtr;
  timestepPtr = _timestepPtr;
  lineNumPtr = _lineNumPtr;
}

void HeaderReader::readHeader() {
  string junk;
  getline(inputStreamPtr->stream, junk);
  inputStreamPtr->stream >> timestepPtr->time;
  // TODO: Read box dimensions for periodic BCs
  for(int i=0; i<7; i++)
    getline(inputStreamPtr->stream, junk);

  *lineNumPtr += 9;

  cout << endl;
  cout << "Timestep " << timestepPtr->time << " @ line " << *lineNumPtr << endl;
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
  int atomId, atomType;
  int ix, iy, iz;

  // Read the line and store data in atom
  // TODO: Account for periodic BCs
  inputStreamPtr->stream >> atomId >> atomType
                         >> atom.x >> atom.y >> atom.z
                         >> ix >> iy >> iz;
  *lineNumPtr++;
  *atomNumPtr++;
}

TimestepReader::TimestepReader() {
  atomNum = 0;
  lineNum = 1;
}
void TimestepReader::setContext(InputStream* _inputStreamPtr, AtomArray* atomArrayPtr, Timestep* _timestepPtr, SimData* _simDataPtr) {
  inputStreamPtr = _inputStreamPtr;
  lineReader.setContext(inputStreamPtr, &atomNum, &lineNum);
  headerReader.setContext(inputStreamPtr, timestepPtr, &lineNum);
  timestepPtr = _timestepPtr;
  simDataPtr = _simDataPtr;
}

void TimestepReader::resetAtomCounter() {
  atomNum = 0;
}

void TimestepReader::readTimestep() {
  resetAtomCounter();
  headerReader.readHeader();
  for(int i=0; i<atomArrayPtr->numAtoms; i++) {
    lineReader.readLine();
    atomArrayPtr->setAtom(i, lineReader.atom);
  }
  timestepPtr->stepNum++;
}

void FrameReader::openStream() {
  inputStream.open(options.dumpfile);
}

void FrameReader::setContext(Options _options, AtomArray* _atomArrayPtr, SimData* _simDataPtr) {
  options = _options;
  atomArrayPtr = _atomArrayPtr;
  simDataPtr = _simDataPtr;

  stepsPerFrame = simDataPtr->stepsPerFrame;
  simDataPtr->setOptions(options);
  openStream();

  timestepReader.setContext(&inputStream, atomArrayPtr, &timestep, simDataPtr);
}

FrameReader::FrameReader() {}

FrameReader::FrameReader(Options options, AtomArray* _atomArrayPtr, SimData* _simDataPtr) {
  setContext(options, _atomArrayPtr, _simDataPtr);
}

void FrameReader::updateFrame() {
  frame.frameStep = timestep.stepNum;
  frame.time = timestep.time;
}

void FrameReader::readFrame() {
  // Read first timestep separately to set frame variables
  // to those of the first timestep in the frame.
  timestepReader.readTimestep();
  updateFrame();

  // Then read the rest
  for(int n; n<stepsPerFrame-1; n++) {
    timestepReader.readTimestep();
  }
}

void InitialTimestepReader::openStream(Options options) {
  inputStream.open(initLoc);
}

void InitialTimestepReader::setContext(Options options, AtomArray* _atomArrayPtr, SimData* _simDataPtr) {
  atomArrayPtr = _atomArrayPtr;
  simDataPtr = _simDataPtr;

  initLoc = "../init.txt";

  openStream(options);

  timestepReader.setContext(&inputStream, atomArrayPtr, &emptyTimestep, simDataPtr);
}

InitialTimestepReader::InitialTimestepReader(Options options, AtomArray* _atomArrayPtr, SimData* _simDataPtr) {
  setContext(options, _atomArrayPtr, _simDataPtr);
}

void InitialTimestepReader::readFromFile() {
  for(int i=0; i<simDataPtr->numAtoms; i++) {
    inputStream.stream >> atom.x >> atom.y >> atom.z;
    atomArrayPtr->setAtom(i, atom);
  }
}

void InitialTimestepReader::readFromStream() {
  // Read first timestep from input file
  // (same file other timesteps will be read from)
  // This will be called if the initial timestep file
  // doesn't exist.
  timestepReader.readTimestep();
}

bool InitialTimestepReader::fileExists() {
  // Return true if file exists, false otherwise.
  ifstream testStream(initLoc);
  return testStream.is_open();
}

void InitialTimestepReader::writeToFile() {
  FILE* outFile = fopen(initLoc.data(), "w");

  for(int i=0; i<simDataPtr->numAtoms; i++) {
    atomArrayPtr->getAtom(i, atom);
    fprintf(outFile, "%10.5f %10.5f %10.5f", atom.x, atom.y, atom.z);
    fflush(outFile);
  }
}

void InitialTimestepReader::readInitialTimestep() {
  // Read first timestep if the file `initLoc` doesn't exist,
  // then write to that file.
  // If it does exist, read from the file
  if(fileExists()) {
    readFromFile();
  }

  else {
    readFromStream();
    writeToFile();
  }
}


///////////////////////
// Top-level Readers //
///////////////////////

DatafileReader::DatafileReader(Options _options, SimData &simData) {
  options = _options;
  simDataPtr = &simData;
  inputStream.open(options.datafile);
  streamPtr = &inputStream.stream;
  read();
}

void DatafileReader::readNumAtoms() {
  bool found;
  int pos = streamPtr->tellg();
  streamPtr->seekg(0, ios::beg);

  found = inputStream.search("atoms");
  *streamPtr >> numAtoms;

  streamPtr->seekg(pos);
  simDataPtr->numAtoms = numAtoms;
  cout << "Read " << numAtoms << " atoms." << endl;
}

void DatafileReader::readMasses() {
  string word;
  int type;
  double mass;
  map<int, bool> foundOne;

  while(!inputStream.nextLineBlank() && streamPtr->good()) {
    *streamPtr >> type >> mass;
    simDataPtr->masses[type] = mass;
  }
}

void DatafileReader::readWaterBonds() {
  int bondId, bondType, OId, HId;
  map<int, bool> foundOne;

  while(!inputStream.nextLineBlank() && streamPtr->good()) {
    // Assume that O atom comes first in bond, then H
    *streamPtr >> bondId >> bondType >> OId >> HId;
    if(bondType == options.waterBondType) {
      if(foundOne[OId]) {
        simDataPtr->waterBonds[OId][1] = HId;
      }
      else {
        simDataPtr->waterBonds[OId][0] = HId;
      }
    }
  }
}

void DatafileReader::read() {
  string section;
  vector<string> sections;

  sections.push_back("Masses");
  sections.push_back("Bonds");

  readNumAtoms();

  for(int i=0; i<2; i++) {
    section = inputStream.searchLine(sections);
    inputStream.skipWhitespace();

    cout << "Reading '" << section  << "'" << endl;
    if(section == "Masses") {
      readMasses();
    }
    else if(section == "Bonds") {
      readWaterBonds();
    }
  }
}

DumpfileReader::DumpfileReader(Options _options, SimData &simData) {
  options = _options;
  simDataPtr = &simData;
  atomArray.setSimData(simData);
  frameReader.setContext(options, &atomArray, &simData);
  countSteps();
}

//Count the number of water atoms in the first timestep
void DumpfileReader::countAtoms() {
  //Count number of atoms
  bool countFlag=true;
  string line;
  int numAtoms=0;
  ifstream *streamPtr;
  int pos;

  streamPtr = &frameReader.inputStream.stream;
  pos = streamPtr->tellg();

  //Ignore the first 9 lines
  for(int i=0;i<9;i++) {
    streamPtr->ignore(256,'\n');
  }

  while(countFlag)
    {
      getline(*streamPtr,line);

      //Count until reaching a line containing "TIMESTEP"
      if(isIn("TIMESTEP", line) || streamPtr->eof()) {
          countFlag=false;
        }
      else
        numAtoms++;
    }

  //Unset eof flag (if set) & return to beginning of file
  streamPtr->clear();
  streamPtr->seekg(pos);

  cout << "Counted " << numAtoms << " atoms." << endl;
  simDataPtr->numAtoms = numAtoms;
}

//Count the number of timesteps
void DumpfileReader::countSteps()
{
  ifstream *streamPtr;
  string line;
  int numSteps=0;
  int lineNum=0;
  int pos;

  streamPtr = &frameReader.inputStream.stream;
  pos = streamPtr->tellg();

  //Count number of timesteps
  while(getline(*streamPtr,line))
    {
      if(isIn("TIMESTEP", line)) {
        numSteps++;
      }
    }

  //Unset eof flag (if set) & return to beginning of file
  streamPtr->clear();
  streamPtr->seekg(pos);

  cout << "Counted " << numSteps << " timesteps." << endl;
  simDataPtr->numSteps = numSteps;
}

void DumpfileReader::readFrame() {
  frameReader.readFrame();
}

bool DumpfileReader::good() {
  return frameReader.frame.frameNum < simDataPtr->numFrames;
}
