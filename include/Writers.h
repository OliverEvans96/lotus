
#ifndef WRITERS_H
#define WRITERS_H

#include "Utils.h"
#include "MDBase.h"
#include "Parameters.h"

using namespace std;

/////////////
// Writers //
/////////////

class FrameWriter
{
  char outputDir[256];
  char path[256];
  FILE* file;
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
  FrameWriter(SimData &simData);
  ~FrameWriter();
  void setOutputDir(char* path);
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

    cout << "QN = " << quantityName << endl;
    cout << "final_fmt = '" << final_fmt << "'" << endl;

    quantityNameArray.push_back(quantityName);
    dataPtrArray.push_back(dataPtr);
    storeType(dataPtr);
    fmtArray.push_back(final_fmt);
    cout << "fmt_arr[" << fmtArray.size()-1 << "] = '" << fmtArray[fmtArray.size()-1] << "'" << endl;
    numQuantities++;
  }

  void deleteFmtStrs();
  void getQuantityStr(char* quantityStr, int i);
  void concatenateQuantityStrs();
  void writeHeader();
  void writeFrame();
};

#endif
