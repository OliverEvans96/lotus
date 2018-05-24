
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
  void openStreams();
};

#endif
