
#ifndef WRITERS_H
#define WRITERS_H

#include "Utils.h"

using namespace std;

/////////////
// Writers //
/////////////

struct FrameWriter
{
  FILE* avgStepData, instStepData;
  void openStreams()
};

#endif
