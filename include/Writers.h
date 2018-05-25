
#ifndef WRITERS_H
#define WRITERS_H

#include "Utils.h"

using namespace std;

/////////////
// Writers //
/////////////

struct FrameWriter
{
  string outputDir;
  map<string, FILE*> streams;
  map<string, string> fmts;
  FILE* avgStepData, instStepData;

  void setOutputDir(string path);
  void openStream(string streamName, string fmt);
  void closeStream(string streamName);
  template <typename T>
  void write(string streamName, T data);
  template <typename T>
  void write(string streamName, vector<T> data, string fmt);
};

#endif
