/**
   @file Readers.h

   Readers for parsing LAMMPS dumpfiles and datafiles.
*/

#ifndef READERS_H
#define READERS_H

#include <cstdio>
#include <map>
#include <cctype>

#include "Utils.h"
#include "MDBase.h"
#include "Atoms.h"
#include "Time.h"

using namespace std;

//////////////////
// Input Stream //
//////////////////

/**
   Convenience wrapper around ifstream.
*/
struct InputStream {
  /// Line number (starting from 1).
  unsigned long int lineNum;
  /// Byte offset in file.
  unsigned long int pos;
  /// Path to file.
  string filename;
  /// Actual ifstream object.
  ifstream stream;

  InputStream();
  /// Open the stream upon instantiation.
  InputStream(string _filename);
  /// Close the stream.
  ~InputStream();
  /// Open the stream and call #verifyStream.
  void open(string _filename);
  /// Check whether #stream is ``good()``.
  void verifyStream();
  /// Skip @p numLines lines forward in the file.
  void skipLines(int numLines);
  /// Consume consecutive whitespace following the cursor.
  void skipWhitespace();
  /// Search for line containing @p term.
  bool search(string term);
  /// Search for line containing an element of @p terms.
  string search(vector<string> terms);
  /// Skip past line containing @p term.
  bool searchLine(string term);
  /// Skip past line containing an element of @p terms.
  string searchLine(vector<string> terms);
  /// Whether the next line is blank.
  bool nextLineBlank();
  /// Read the next line and return the cursor to its previous position.
  string peekLine();
};

/////////////
// Readers //
/////////////

/**
   Reads the header from a LAMMPS dumpfile.
   Read the timestep and box bounds.
   @note Could also (but currently does not) read the number of atoms.

   @see TimestepReader.
*/
class HeaderReader {
  Options options;
  InputStream* inputStreamPtr;
  Timestep* timestepPtr;
  BoxBounds* boundsPtr;
  SimData* simDataPtr;

 public:
  int* lineNumPtr;
  /** Copy #options and save pointers to relevant quantities.
      Called by TimestepReader::setContext.
  */
  void setContext(Options options, InputStream* _inputStreamPtr, SimData* _simDataPtr, Timestep* _timestepPtr, int* _lineNumPtr, BoxBounds* _boundsPtr);
  void readHeader();
};

/**
   Reads one line at a time of atom positions from LAMMPS dumpfile.
   Reads scaled positions and calculates unitful positions.

   @see TimestepReader.
*/
class LineReader {
  Options options;
  string line;
  string input;
  int* atomNumPtr;
  int* lineNumPtr;
  InputStream* inputStreamPtr;
  BoxBounds* boundsPtr;

 public:
  Atom atom;
  /** Copy #options and save pointers to relevant quantities.
      Called by TimestepReader::setContext.
  */
  void setContext(Options options, InputStream *_inputStreamPtr, int* _atomNumPtr, int* _lineNumPtr, BoxBounds* _boundsPtr);
  void readLine();
};

/**
   Read one timestep at a time from a LAMMP dumpfile.
   Reads both the header and position lines.
   @see FrameReader.
*/
class TimestepReader {
  Options options;
  LineReader lineReader;
  HeaderReader headerReader;
  Timestep* timestepPtr;
  AtomArray* atomArrayPtr;
  InputStream* inputStreamPtr;
  SimData* simDataPtr;
  BoxBounds timestepBounds;
  double xlo, xhi;
  double ylo, yhi;
  double zlo, zhi;

 public:
  int atomNum;
  int lineNum;

  TimestepReader();
  /** Copy #options and save pointers to relevant quantities.
      Called by FrameReader::setContext.
  */
  void setContext(Options _options, InputStream* _inputStreamPtr, AtomArray* _atomArrayPtr, Timestep* _timestepPtr, SimData* _simDataPtr);
  /// Set #atomNum = 0.
  void resetAtomCounter();
  /// Read one timestep.
  void readTimestep(int stepInFrame);
};

/**
   Read one frame, containing multiple timesteps from a LAMMP dumpfile.
   The size of a frame is determined by Options::stepsPerFrame
   (although the last frame may be a different size).
   @see Frame::stepsThisFrame.
   @see LastFrame::numSteps.
   @see DumpfileReader.
*/
class FrameReader {
 public:

  Options options;
  TimestepReader timestepReader;
  AtomArray* atomArrayPtr;
  Timestep timestep;
  SimData* simDataPtr;
  InputStream inputStream;
  Frame frame;
  int stepsPerFrame;

  FrameReader();
  FrameReader(Options options, AtomArray* _atomArrayPtr, SimData* _simDataPtr);
  /** Copy #options and save pointers to relevant quantities.
      Called by DumpfileReader::DumpfileReader.
  */
  void setContext(Options options, AtomArray* _atomArrayPtr, SimData* _simDataPtr);
  /// Open #inputStream with Options::dumpfile.
  void openStream();
  /// Set frame variables from first timestep in frame.
  void updateFrame();
  /// Read one frame.
  void readFrame();

  /// Count the number of atoms in the simulation (in one timestep).
  void countAtoms();
  /** Count the number of timesteps in the dumpfile.
      Called by DumpfileReader::DumpfileReader.
  */
  void countSteps();
};

/**
   Read the first timestep in a simulation (not just one dumpfile).
   Read from dumpfile or from another dedicated file
   if this is not the first dumpfile.

   This is useful for calculating dynamical quantities such as MSD.

   This is not currently used anywhere.
*/
class InitialTimestepReader {
  Options options;
 public:
  InitialTimestepReader(Options options, AtomArray* _atomArrayPtr, SimData* _simDataPtr);
  void setContext(Options options, AtomArray* _atomArrayPtr, SimData* _simDataPtr);
  string initLoc;
  TimestepReader timestepReader;
  AtomArray* atomArrayPtr;
  Timestep* timestepPtr;
  InputStream inputStream;
  Timestep emptyTimestep;
  SimData* simDataPtr;
  Atom atom;

  bool fileExists();
  void openStream(Options _options);
  void readFromFile();
  void readFromStream();
  void writeToFile();
  void readInitialTimestep();
};

///////////////////////
// Top-level Readers //
///////////////////////

/**
   Read a LAMMPS Datafile.

   Reads box dimensions, number of atoms, masses of each type,
   and water bonds (useful if water molecule orientation is desired).
   @note Water bonds are not currently utilized.
*/
class DatafileReader {
  Options options;
  InputStream inputStream;
  vector<int> liquidTypes;
  /// Atom type of water hydrogen atoms
  int HType;
  /// Atom type of water oxygen atoms
  int OType;
  ifstream *streamPtr;
  SimData *simDataPtr;
  int numAtoms;

  void readNumAtoms();
  void readMasses();
  void readWaterBonds();
  void readBoxBounds();

 public:
  DatafileReader(SimData &simData);
  /// Open Options::datafile for reading.
  void openStream();
  /// Read the whole datafile.
  void read();
};

/**
   Read a LAMMPS Dumpfile one frame at a time.
*/
struct DumpfileReader {
  Options options;
  FrameReader frameReader;
  SimData *simDataPtr;
  AtomArray *atomArrayPtr;
  /** Index of current frame.
     @see Frame::frameNum.
  */
  int frameNum;

  DumpfileReader(AtomArray &atomArray);
  /// Count the number of water atoms in the first timestep
  void countAtoms();
  void countSteps();
  void readFrame();
  /// Whether another frame can be read.
  bool good();
};


#endif
