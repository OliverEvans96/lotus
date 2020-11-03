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
  stream.open(filename.data());
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

/**
   Skip whitespace in input stream, including newlines,
   until the next non-whitespace character.
**/
void InputStream::skipWhitespace() {
  int num;
  while(isspace(stream.peek())) {
    stream.get();
  }
}

/**
   Move the ifstream cursor to the beginning of
   the first line containing @p term.
   (the matching line will not be consumed)
   Returns whether the line was found
*/
bool InputStream::search(string term) {
  string line;
  bool found;
  int pos;

  found = false;

  while(stream.good()) {
    // Save beginning of line
    pos = stream.tellg();
    getline(stream, line);
    lineNum++;
    if(isIn(term, line)) {
      found = true;
      // Rewind to line beginning
      stream.seekg(pos);
      break;
    }
  }

  return found;
}

/**
   Move the ifstream cursor one line past
   the first line containing one of the
   strings in @p terms.
   (the matching line will be consumed)
   Returns the line found
*/
string InputStream::search(vector<string> terms) {
  string line;
  string term;
  bool found;
  vector<string>::iterator it;

  term = "";
  found = false;

  while(stream.good()) {
    getline(stream, line);
    lineNum++;
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

/**
   Move the ifstream cursor one line past
   the first line containing @p term.
   (the matching line will be consumed)
   Returns whether the line was found.
*/
bool InputStream::searchLine(string term) {

  string line;
  bool found;

  found = false;

  while(stream.good()) {
    getline(stream, line);
    lineNum++;
    if(isIn(term, line)) {
      found = true;
      break;
    }
  }

  if(!found) {
    cout << "SearchLine(string) failed for term '" << term << "'" << endl;
  }

  return found;
}

/*
   Move the ifstream cursor one line past
   the first line containing one of the
   strings in `terms`.
   (the matching line will be consumed)
   Returns the term found
*/
string InputStream::searchLine(vector<string> terms) {

  string line;
  string term;
  bool found;
  vector<string>::iterator it;

  found = false;

  while(!found && stream.good()) {
    getline(stream, line);
    lineNum++;
    for(it=terms.begin(); it<terms.end(); it++) {
      term = *it;
      if(isIn(term, line)) {
        found = true;
        break;
      }
    }
  }

  if(!found) {
    cout << "SearchLine(vector<string>) failed for terms: ";
    for(it=terms.begin(); it<terms.end(); it++) {
      term = *it;
      cout << "'" << term << "' ";
    }
    cout << endl;
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
  pos = stream.tellg();
  getline(stream, line);
  stream.seekg(pos);
  return line;
}

/////////////
// Readers //
/////////////

void HeaderReader::setContext(Options _options, InputStream* _inputStreamPtr, SimData* _simDataPtr, Timestep* _timestepPtr, int* _lineNumPtr, BoxBounds* _boundsPtr) {
  options = _options;
  inputStreamPtr = _inputStreamPtr;
  simDataPtr = _simDataPtr;
  timestepPtr = _timestepPtr;
  lineNumPtr = _lineNumPtr;
  boundsPtr = _boundsPtr;
}

void HeaderReader::readHeader() {
  string junk;
  if(options.verbose)
    cout << "Reading header: @ " << inputStreamPtr->stream.tellg() << " '" << inputStreamPtr->peekLine() << "'" << endl;

  getline(inputStreamPtr->stream, junk); // ITEM: TIMESTEP
  inputStreamPtr->stream >> timestepPtr->time;
  inputStreamPtr->stream.ignore(256, '\n');

  getline(inputStreamPtr->stream, junk); // ITEM: NUMBER OF ATOMS
  inputStreamPtr->stream >> simDataPtr->numAtoms;
  inputStreamPtr->stream.ignore(256, '\n');

  getline(inputStreamPtr->stream, junk); // ITEM: BOX BOUNDS pp pp pp

  inputStreamPtr->stream >> boundsPtr->xlo >> boundsPtr->xhi;
  inputStreamPtr->stream >> boundsPtr->ylo >> boundsPtr->yhi;
  inputStreamPtr->stream >> boundsPtr->zlo >> boundsPtr->zhi;
  inputStreamPtr->stream.ignore(256, '\n');

  getline(inputStreamPtr->stream, junk); // ITEM: ATOMS type x y z

  *(lineNumPtr) += 9;

  cout << endl;
  cout << "Timestep " << timestepPtr->time << " @ line " << *lineNumPtr << endl;
}

void LineReader::setContext(Options _options, InputStream* _inputStreamPtr, int* _atomNumPtr, int* _lineNumPtr, BoxBounds* _boundsPtr) {
  options = _options;
  inputStreamPtr = _inputStreamPtr;
  atomNumPtr = _atomNumPtr;
  lineNumPtr = _lineNumPtr;
  boundsPtr = _boundsPtr;

  // Chose which line reading function to use
  switch(options.lineFormat) {
  case 1:
    _readLineFunctionPtr = &LineReader::readLineFormat1;
    break;
  case 2:
    _readLineFunctionPtr = &LineReader::readLineFormat2;
    break;
  default:
    cout << "ERROR: INVALID LINE FORMAT " << options.lineFormat << endl;
    throw 1;
  }
}

void LineReader::readLine() {
  // Call the function pointed to by the pointer
  (this->*_readLineFunctionPtr)();
}

/// Read the line (format 1) and store data in #atom.
void LineReader::readLineFormat1() {
  int atomId;
  double xs, ys, zs;
  int ix, iy, iz;

  // xs, ys, zs are between 0 and 1,
  // ix, iy, iz are number of periodic image
  inputStreamPtr->stream >> atomId >> atom.type
                         >> xs >> ys >> zs
                         >> ix >> iy >> iz;

  // Ignore ix, iy, iz so that all
  // periodic images are in same box
  atom.x = boundsPtr->xlo + (boundsPtr->xhi - boundsPtr->xlo) * xs - options.xCenter;
  atom.y = boundsPtr->ylo + (boundsPtr->yhi - boundsPtr->ylo) * ys - options.yCenter;
  atom.z = boundsPtr->zlo + (boundsPtr->zhi - boundsPtr->zlo) * zs;

  inputStreamPtr->stream.ignore(256, '\n');
  (*lineNumPtr)++;
  (*atomNumPtr)++;
}


/// Read the line (format 2) and store data in #atom.
void LineReader::readLineFormat2() {
  double x, y;
  inputStreamPtr->stream >> atom.type >> x >> y >> atom.z;
  atom.x = x - options.xCenter;
  atom.y = y - options.yCenter;
  inputStreamPtr->stream.ignore(256, '\n');
  (*lineNumPtr)++;
  (*atomNumPtr)++;
}

TimestepReader::TimestepReader() {
  atomNum = 0;
  lineNum = 1;
}
void TimestepReader::setContext(Options _options, InputStream* _inputStreamPtr, AtomArray* _atomArrayPtr, Timestep* _timestepPtr, SimData* _simDataPtr) {
  options = _options;
  inputStreamPtr = _inputStreamPtr;
  timestepPtr = _timestepPtr;
  simDataPtr = _simDataPtr;
  atomArrayPtr = _atomArrayPtr;
  lineReader.setContext(options, inputStreamPtr, &atomNum, &lineNum, &timestepBounds);
  headerReader.setContext(options, inputStreamPtr, simDataPtr, timestepPtr, &lineNum, &timestepBounds);
}

void TimestepReader::resetAtomCounter() {
  atomNum = 0;
}

/**
   @p stepInFrame is necessary to determine the index
   to store atom position in #atomArrayPtr.
*/
void TimestepReader::readTimestep(int stepInFrame) {
  if(options.verbose)
    cout << "Reading timestep: @ " << inputStreamPtr->stream.tellg() << " '" << inputStreamPtr->peekLine() << "'" << endl;
  resetAtomCounter();
  headerReader.readHeader();
  for(int atomNum=0; atomNum<atomArrayPtr->numAtoms; atomNum++) {
    lineReader.readLine();
    // All atoms in frame are stored,
    // so step offset must be considered.
    atomArrayPtr->setAtom(atomNum, stepInFrame, lineReader.atom);
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

  // Save pointer to frame for simData
  simDataPtr->framePtr = &frame;

  timestepReader.setContext(options, &inputStream, atomArrayPtr, &timestep, simDataPtr);
}

FrameReader::FrameReader() {}

FrameReader::FrameReader(Options options, AtomArray* _atomArrayPtr, SimData* _simDataPtr) {
  setContext(options, _atomArrayPtr, _simDataPtr);
}

/**
   Timestep::stepNum -> Frame::frameStep
   Timestep::time -> Frame::time
*/
void FrameReader::updateFrame() {
  frame.frameStep = timestep.stepNum;
  frame.time = timestep.time;
}

/**
   Read first timestep separately to set frame variables
   to those of the first timestep in the frame.
*/
void FrameReader::readFrame() {
  if(options.verbose)
    cout << "Reading frame: @ " << inputStream.stream.tellg() << " '" << inputStream.peekLine() << "'" << endl;
  timestepReader.readTimestep(0);
  updateFrame();

  if(frame.frameNum < simDataPtr->numFrames-1) {
    frame.stepsThisFrame = stepsPerFrame;
  }
  else {
    frame.stepsThisFrame = simDataPtr->lastFrame.numSteps;
  }

  frame.atomsThisFrame = simDataPtr->numAtoms * frame.stepsThisFrame;

  if(options.verbose) {
    cout << "frameNum = " << frame.frameNum << endl;
    cout << "numFrames = " << simDataPtr->numFrames << endl;
    cout << "stepsThisFrame = " << frame.stepsThisFrame << endl;
  }

  // Then read the rest
  for(int n=1; n<frame.stepsThisFrame; n++) {
    if(options.verbose)
      cout << "frameStep " << n << endl;
    timestepReader.readTimestep(n);
  }

  frame.frameNum++;
}

///////////////////////
// Top-level Readers //
///////////////////////

DatafileReader::DatafileReader(SimData &simData) {
  options = simData.options;
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

void DatafileReader::readBoxBounds() {
  bool found;
  int pos = streamPtr->tellg();
  cout << "bounds initial pos=" << pos << endl;
  streamPtr->seekg(0);
  found = inputStream.search("xlo");
  if(options.verbose) {
    cout << "Found x: " << inputStream.peekLine() << endl;
  }
  *streamPtr >> simDataPtr->simBounds.xlo >> simDataPtr->simBounds.xhi;

  streamPtr->seekg(0);
  found = inputStream.search("ylo");
  if(options.verbose) {
    cout << "Found y: " << inputStream.peekLine() << endl;
  }
  *streamPtr >> simDataPtr->simBounds.ylo >> simDataPtr->simBounds.yhi;

  streamPtr->seekg(0);
  found = inputStream.search("zlo");
  if(options.verbose) {
    cout << "Found z: " << inputStream.peekLine() << endl;
  }
  *streamPtr >> simDataPtr->simBounds.zlo >> simDataPtr->simBounds.zhi;

  cout << "Read bounds:" << endl;
  cout << "x: " << simDataPtr->simBounds.xlo << " " << simDataPtr->simBounds.xhi << endl;
  cout << "y: " << simDataPtr->simBounds.ylo << " " << simDataPtr->simBounds.yhi << endl;
  cout << "z: " << simDataPtr->simBounds.zlo << " " << simDataPtr->simBounds.zhi << endl;

  streamPtr->seekg(pos);
}

void DatafileReader::readMasses() {
  string word;
  int type;
  double mass;
  map<int, bool> foundOne;

  while(!inputStream.nextLineBlank() && streamPtr->good()) {
    *streamPtr >> type >> mass;
    simDataPtr->masses[type] = mass;
    // Ignore the rest of the line (e.g. comments)
    streamPtr->ignore(256,'\n');
    // Back up one character
    streamPtr->seekg(-1,ios_base::cur);
    inputStream.lineNum++;
  }
}

/// Assume that O atom comes first in bond, then H
void DatafileReader::readWaterBonds() {
  int bondId, bondType, OId, HId;
  map<int, bool> foundOne;

  while(!inputStream.nextLineBlank() && streamPtr->good()) {
    *streamPtr >> bondId >> bondType >> OId >> HId;
    if(bondType == options.waterBondType) {
      if(foundOne[OId]) {
        simDataPtr->waterBonds[OId][1] = HId;
      }
      else {
        simDataPtr->waterBonds[OId] = new int[2];
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

  if(simDataPtr->options.numAtoms <= 0) {
    readNumAtoms();
  }
  readBoxBounds();

  // Search through the data file for these two sections
  for(int i=0; i<sections.size(); i++) {
    section = inputStream.searchLine(sections);
    inputStream.skipWhitespace();

    cout << "Reading '" << section << "' @ " << inputStream.stream.tellg() << " '" << inputStream.peekLine() << "'" << endl;
    if(section == "Masses") {
      readMasses();
    }
    else if(section == "Bonds") {
      readWaterBonds();
    }
  }
}

DumpfileReader::DumpfileReader(AtomArray &atomArray) {
  atomArrayPtr = &atomArray;
  simDataPtr = atomArrayPtr->simDataPtr;
  options = simDataPtr->options;
  frameReader.setContext(options, atomArrayPtr, simDataPtr);
  countSteps();
}

/// @note This method is not currently called anywhere.
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

  cout << "Counting timesteps" << endl;

  // 9 lines for header, rest are atoms
  int linesPerStep = simDataPtr->numAtoms + 9;

  streamPtr = &frameReader.inputStream.stream;
  pos = streamPtr->tellg();

  //Count number of timesteps
  while(getline(*streamPtr,line)) {
      lineNum++;
  }

  numSteps = lineNum / linesPerStep;

  //Unset eof flag (if set) & return to beginning of file
  streamPtr->clear();
  streamPtr->seekg(pos);

  simDataPtr->setNumSteps(numSteps);
  atomArrayPtr->allocateArrays();
  cout << "Counted " << numSteps << " timesteps." << endl;
}

void DumpfileReader::readFrame() {
  frameNum = frameReader.frame.frameNum;
  frameReader.readFrame();

  cout << "frameNum = " << frameNum << endl;
  if(options.verbose) {
    atomArrayPtr->printStats();
  }
}

/// Compare frame numbers to see if this is the end.
bool DumpfileReader::good() {
  return frameReader.frame.frameNum < simDataPtr->numFrames;
}
