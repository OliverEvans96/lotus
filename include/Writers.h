
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

class WriterBase {
 protected:
  char outputDir[256];
  char path[256];
  SimData *simDataPtr;
  Options options;
  bool leftJustify;

  // Maximum line length is 2047 characters (last element is null).
  char line[2048];

  void storeType(double *x);
  void storeType(int *x);

  void setOutputDir(const char* path);
  void createDataDir();

  void getFmtStr(char* fmt, double* dataPtr);
  void getFmtStr(char* fmt, int* dataPtr);
  void getFmtStr(char* fmt, const char* dataPtr);
  void deleteFmtStrs();

  int numQuantities;

  vector<void*> dataPtrArray;
  vector<char> typeArray;
  vector<const char*> quantityNameArray;
  vector<const char*> fmtArray;

 public:
  WriterBase(SimData &simData);
  ~WriterBase();
};

class ScalarWriter : public WriterBase {
  vector<FILE*> files;
  DumpfileReader *dumpfileReaderPtr;
  Droplet *dropletPtr;

 public:
  ScalarWriter(DumpfileReader &dumpfileReader, Droplet &droplet);
  ~ScalarWriter();

  FILE* openFileBase(const char* filename);
  void openFile(const char* quantityName);
  void closeFiles();

  void setOutputQuantities();
  void getQuantityStr(char* quantityStr, int i);
  void writeHeaders();
  void writeFrame();

  // For some reason, this template function must
  // be fully defined in header file.
  // http://www.cplusplus.com/forum/general/114299/
  // https://stackoverflow.com/questions/10632251/undefined-reference-to-template-function
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

class ArrayWriter : public WriterBase {
  void setOutputQuantities();
  void getQuantityStr(char* quantityStr, int quantityNum, int i, int j);
  void concatenateQuantityStrs(int quantityNum, int i);
  void createQuantityDir(const char* quantityName);
  void getFilePath(char* filePath, const char* quantityName);
  void writeHeader(FILE* file, int quantityNum);
  void writeQuantityFile(int i);

  vector<int> numColumnsArray;
  vector<int*> lengthPtrArray;
  vector<char**> headersArray;
  DumpfileReader *dumpfileReaderPtr;
  Droplet *dropletPtr;

 public:
  ArrayWriter(DumpfileReader &dumpfileReader, Droplet &droplet);
  ~ArrayWriter();
  void writeHeaders();
  void writeFrame();

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
