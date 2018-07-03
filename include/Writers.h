/**
   @file Writers.h

   Classes for writing time-dependent text output files.

   @rst
   These classes produce the files in the ``data`` directory
   For a description of the output files, see :ref:`output.rst`.
   @endrst
*/

#ifndef WRITERS_H
#define WRITERS_H

#include "Utils.h"
#include "MDBase.h"
#include "Parameters.h"
#include "Droplet.h"
#include "Substrate.h"
#include "Readers.h"

using namespace std;

/////////////
// Writers //
/////////////

/**
   Common functions for scalar and array writers.

   This class is not instantiated directly, only it's children are.

   @see ScalarWriter
   @see ArrayWriter.
*/
class WriterBase {
 protected:
  /// @see Options::outLoc
  char outputDir[256];
  char path[256];
  SimData *simDataPtr;
  Options options;

  /// Whether to left-justify output data
  bool leftJustify;

  /// Number of quantities to be written
  int numQuantities;

  /// Pointers to actual output quantity data
  vector<void*> dataPtrArray;
  /// `char`s representing output quantity types
  vector<char> typeArray;
  /// Labels for output quantities
  vector<const char*> quantityNameArray;
  /// Format strings (e.g. '%15.6f') for output quantities
  vector<const char*> fmtArray;

  /// Maximum line length is 2047 characters (last element is null).
  char line[2048];

  /// Store a char (d) representing the type of @p *x in @p typeArray.
  void storeType(double *x);
  /// Store a char (i) representing the type of @p *x in @p typeArray.
  void storeType(int *x);

  void setOutputDir(const char* path);
  /// Create data directory if it doesn't already exist
  void createDataDir();

  /// Get format string, e.g. '%15.6f'.
  void getFmtStr(char* fmt, double* dataPtr);
  /// Get format string, e.g. '%15s'.
  void getFmtStr(char* fmt, int* dataPtr);
  /// Get format string, e.g. '%15s'.
  void getFmtStr(char* fmt, const char* dataPtr);
  void deleteFmtStrs();

 public:
  WriterBase(SimData &simData);
  ~WriterBase();
};

/**
   Write single column scalar data.

   @rst
   See :ref:`output-data`.
   @endrst
*/
class ScalarWriter : public WriterBase {
  /// All of the files being written to each timestep.
  vector<FILE*> files;
  DumpfileReader *dumpfileReaderPtr;
  Droplet *dropletPtr;

 public:
  ScalarWriter(DumpfileReader &dumpfileReader, Droplet &droplet);
  ~ScalarWriter();

  /// Open a file at Options::outLoc/data/`filename`.
  FILE* openFileBase(const char* filename);
  /// Call #openFileBase and store the @p FILE* in #files
  void openFile(const char* quantityName);
  /// Close all `FILE*`s in #files.
  void closeFiles();

  /** Call #addQuantity for each quantity of interest.
      If you'd like to add a new output file
      or modify file names,
      this is the place to do so.
  */
  void setOutputQuantities();
  /// Get the formatted data for this quantity.
  void getQuantityStr(char* quantityStr, int i);
  /** Write the quantity name to the first line of each file,
      with '#' prepended.
  */
  void writeHeaders();
  /// Write the values of all scalar quantities to the appropriate files.
  void writeFrame();

  /**
     Register a scalar quantity to be written to an output file.

     A pointer to the data must be given since the value
     will presumably change over time.

     @param quantityName   Name of the quantity
     @param dataPtr        Pointer to the data

     @note For some reason, this template function must
     be fully defined in header file.

     @see http://www.cplusplus.com/forum/general/114299/
     @see https://stackoverflow.com/questions/10632251/undefined-reference-to-template-function
  */
  template <typename T>
  void addQuantity(const char* quantityName, T *dataPtr) {
    char* fmt = new char[16];
    quantityNameArray.push_back(quantityName);
    openFile(quantityName);
    dataPtrArray.push_back(dataPtr);
    storeType(dataPtr);
    getFmtStr(fmt, dataPtr);
    fmtArray.push_back(fmt);
    numQuantities++;
  }
};

/**
   Write multi-column array data.

   @rst
   See :ref:`output-data`.
   @endrst
*/
class ArrayWriter : public WriterBase {
  DumpfileReader *dumpfileReaderPtr;
  Droplet *dropletPtr;
  /// The number of columns in each quantity.
  vector<int> numColumnsArray;
  /** A pointer to the number of rows of data,
      which may change from one step to the next.
  */
  vector<int*> lengthPtrArray;
  /** A pointer to an array of headers whose length is defined
      by the corresponding element of #numColumnsArray.
  */
  vector<char**> headersArray;
  /** Call #addQuantity for each quantity of interest.
      If you'd like to add a new output file
      or modify file names,
      this is the place to do so.
  */

  void setOutputQuantities();
  /// Get the formatted data for one element in the array quantity.
  void getQuantityStr(char* quantityStr, int quantityNum, int i, int j);
  /// Join individual quantity strings to form a row.
  void concatenateQuantityStrs(int quantityNum, int i);
  /// Create Options::outLoc/data/`quantityName` if it doesn't alread exist.
  void createQuantityDir(const char* quantityName);
  /// Get the path to the quantity file for the current timestep.
  void getFilePath(char* filePath, const char* quantityName);
  /// Write column headers to the first line of the file.
  void writeHeader(FILE* file, int quantityNum);
  /** Write the array data for quantity @p i
      to the appropriate file for this timestep.
  */
  void writeQuantityFile(int i);

 public:
  ArrayWriter(DumpfileReader &dumpfileReader, Droplet &droplet);
  ~ArrayWriter();
  void writeHeaders();
  void writeFrame();

  /**
     Register an array quantity to be written
     to a new output file each timestep.

     Pointers to the data and number of rows must be given
     since their values will presumably change over time.

     @param quantityName   Name of the quantity.
     @param dataPtr        Pointer to the data
                           (may vary over time).
     @param lengthPtr      Pointer to the number of rows in the array
                           (may vary over time).
     @param headers        Pointer to array of column headers
                           (e.g. {"x", "y"}).
     @param numColumns     Number of columns in the array
                           (assumed constant over time).

     @note For some reason, this template function must
     be fully defined in header file.

     @see http://www.cplusplus.com/forum/general/114299/
     @see https://stackoverflow.com/questions/10632251/undefined-reference-to-template-function
  */
  template <typename T>
  void addQuantity(const char* quantityName, T **dataPtr, int *lengthPtr, char** headers, int numColumns) {
    char* fmt = new char[16];

    quantityNameArray.push_back(quantityName);
    createQuantityDir(quantityName);
    dataPtrArray.push_back(dataPtr);
    lengthPtrArray.push_back(lengthPtr);
    headersArray.push_back(headers);
    numColumnsArray.push_back(numColumns);
    storeType(*dataPtr);
    getFmtStr(fmt, *dataPtr);
    fmtArray.push_back(fmt);
    numQuantities++;
  }
};

#endif
