/**
   @file Readers.h

   @brief Readers for parsing LAMMPS dumpfiles and datafiles.
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
   @brief Convenience wrapper around ifstream.
*/
struct InputStream {
  /// @brief Line number (starting from 1).
  unsigned long int lineNum;
  /// @brief Byte offset in file.
  unsigned long int pos;
  /// @brief Path to file.
  string filename;
  /// @brief Actual ifstream object.
  ifstream stream;

  InputStream();
  /// @brief Open the stream upon instantiation.
  InputStream(string _filename);
  /// @brief Close the stream.
  ~InputStream();
  /// @brief Open the stream and call #verifyStream.
  void open(string _filename);
  /// @brief Check whether #stream is ``good()``.
  void verifyStream();
  /// @brief Skip @p numLines lines forward in the file.
  void skipLines(int numLines);
  /// @brief Consume consecutive whitespace following the cursor.
  void skipWhitespace();
  /// @brief Search for line containing @p term.
  bool search(string term);
  /// @brief Search for line containing an element of @p terms.
  string search(vector<string> terms);
  /// @brief Skip past line containing @p term.
  bool searchLine(string term);
  /// @brief Skip past line containing an element of @p terms.
  string searchLine(vector<string> terms);
  /// @brief Whether the next line is blank.
  bool nextLineBlank();
  /// @brief Read the next line and return the cursor to its previous position.
  string peekLine();
};

/////////////
// Readers //
/////////////

/**
   @brief Reads the header from a LAMMPS dumpfile.

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
  /** @brief Copy #options and save pointers to relevant quantities.

      Called by TimestepReader::setContext.
  */
  void setContext(Options options, InputStream* _inputStreamPtr, SimData* _simDataPtr, Timestep* _timestepPtr, int* _lineNumPtr, BoxBounds* _boundsPtr);
  void readHeader();
};

/**
   @brief Reads one line at a time of atom positions from LAMMPS dumpfile.

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
  /** @brief Copy #options and save pointers to relevant quantities.

      Called by TimestepReader::setContext.
  */
  void setContext(Options options, InputStream *_inputStreamPtr, int* _atomNumPtr, int* _lineNumPtr, BoxBounds* _boundsPtr);
  void readLine();
};

/**
   @brief Read one timestep at a time from a LAMMP dumpfile.

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
  /** @brief Copy #options and save pointers to relevant quantities.

      Called by FrameReader::setContext.
  */
  void setContext(Options _options, InputStream* _inputStreamPtr, AtomArray* _atomArrayPtr, Timestep* _timestepPtr, SimData* _simDataPtr);
  /// @brief Set #atomNum = 0.
  void resetAtomCounter();
  /// @brief Read one timestep.
  void readTimestep(int stepInFrame);
};

/**
   @brief Read one frame, containing multiple timesteps from a LAMMP dumpfile.

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
  /** @brief Copy #options and save pointers to relevant quantities.

      Called by DumpfileReader::DumpfileReader.
  */
  void setContext(Options options, AtomArray* _atomArrayPtr, SimData* _simDataPtr);
  /// @brief Open #inputStream with Options::dumpfile.
  void openStream();
  /// @brief Set frame variables from first timestep in frame.
  void updateFrame();
  /// @brief Read one frame.
  void readFrame();

  /// @brief Count the number of atoms in the simulation (in one timestep).
  void countAtoms();
  /** @brief Count the number of timesteps in the dumpfile.

      Called by DumpfileReader::DumpfileReader.
  */
  void countSteps();
};

///////////////////////
// Top-level Readers //
///////////////////////

/**
   @brief Read a LAMMPS Datafile.

   Reads box dimensions, number of atoms, masses of each type,
   and water bonds (useful if water molecule orientation is desired).
   @note Water bonds are not currently utilized.
*/
class DatafileReader {
  Options options;
  InputStream inputStream;
  vector<int> liquidTypes;
  /// @brief Atom type of water hydrogen atoms
  int HType;
  /// @brief Atom type of water oxygen atoms
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
  /// @brief Open Options::datafile for reading.
  void openStream();
  /// @brief Read the whole datafile.
  void read();
};

/**
   @brief Read a LAMMPS Dumpfile one frame at a time.
*/
struct DumpfileReader {
  Options options;
  FrameReader frameReader;
  SimData *simDataPtr;
  AtomArray *atomArrayPtr;
  /** @brief Index of current frame.
     @see Frame::frameNum.
  */
  int frameNum;

  DumpfileReader(AtomArray &atomArray);
  /// @brief Count the number of water atoms in the first timestep
  void countAtoms();
  void countSteps();
  void readFrame();
  /// @brief Whether another frame can be read.
  bool good();
};

#endif
