
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

class FrameWriter
{
  char outputDir[256];
  char path[256];
  FILE* file;
  DumpfileReader *dumpfileReaderPtr;
  Droplet *dropletPtr;
  SimData *simDataPtr;
  Options options;

  vector<const char*> quantityNameArray;
  vector<void*> dataPtrArray;
  vector<char> typeArray;
  vector<const char*> fmtArray;
  int numQuantities = 0;

  // Maximum line length is 2047 characters (last element is null).
  char line[2048];

  void storeType(double *x);
  void storeType(int *x);

  void stripTrailingSlash(char* strippedPath, char* path);
  void joinPath(char* path, char* prefix, char* suffix);

 public:
  FrameWriter(char* filename, DumpfileReader &dumpfileReader, Droplet &droplet);
  ~FrameWriter();

  void setOutputDir(char* path);
  void setOutputQuantities();
  void openFile(char* filename);
  void closeFile();

  void getDefaultFmt(char* final_fmt, double* dataPtr);
  void getDefaultFmt(char* final_fmt, int* dataPtr);

  // For some reason, this template function must
  // be fully defined in header file.
  // http://www.cplusplus.com/forum/general/114299/
  // https://stackoverflow.com/questions/10632251/undefined-reference-to-template-function
  template <typename T>
  void addQuantity(const char* quantityName, T *dataPtr, const char* fmt=" ") {
    char* final_fmt = new char[16];
    // Get default format string based on variable type
    if(strcmp(fmt, " ") == 0) {
      getDefaultFmt(final_fmt, dataPtr);
    }
    else {
      strcpy(final_fmt, fmt);
    }

    quantityNameArray.push_back(quantityName);
    dataPtrArray.push_back(dataPtr);
    storeType(dataPtr);
    fmtArray.push_back(final_fmt);
    numQuantities++;
  }

  void deleteFmtStrs();
  void getQuantityStr(char* quantityStr, int i);
  void concatenateQuantityStrs();
  void writeHeader();
  void writeFrame();
};

#endif
