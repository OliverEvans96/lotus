
#include "Writers.h"

/////////////
// Writers //
/////////////

FrameWriter::openStreams() {
  //Output file streams
  FILE* avgStepData = fopen("avgStepData.txt","w");
  fprintf(avgStepData,"%8.8s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s\n","#time","bulkEdge","monoEdge","contactAngle","dropletHeight","avgRadius","avgDipole","MSDx","MSDy","MSDz","MSDavg","frameFlux","dMe_dt","nMono","chi2s");
  fflush(avgStepData);

  FILE* instStepData = fopen("instStepData.txt","w");
  fprintf(instStepData,"%8.8s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s \n","#time","avgRadius","avgDipole","MSDx","MSDy","MSDz","MSDavg","stepFlux");
  fflush(instStepData);
}

FrameWriter::FrameWriter() {
  openStreams();
}
