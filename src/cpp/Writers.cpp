
#include "Writers.h"

/////////////
// Writers //
/////////////

FrameWriter::FrameWriter(const char* filename, DumpfileReader &dumpfileReader, Droplet &droplet) {
  dumpfileReaderPtr = &dumpfileReader;
  dropletPtr = &droplet;
  simDataPtr = dumpfileReaderPtr->simDataPtr;
  options = simDataPtr->options;

  setOutputDir(options.outLoc.data());
  openFile(filename);

  setOutputQuantities();
  writeHeader();
}
FrameWriter::~FrameWriter() {
  closeFile();
  deleteFmtStrs();
}

void FrameWriter::storeType(double *x) {
  typeArray.push_back('d');
}
void FrameWriter::storeType(int *x) {
  typeArray.push_back('i');
}

void FrameWriter::setOutputDir(const char* path) {
  strcpy(outputDir, path);
  if(options.verbose) {
    cout << "Set output dir '" << outputDir << "'" << endl;
  }
}

// Choose which quantities are written to the main results file
// A quantity name and pointer to the quantity are given for each
void FrameWriter::setOutputQuantities() {
  addQuantity("timestep", &dumpfileReaderPtr->frameReader.frame.time);
  addQuantity("monoEdge", &dropletPtr->monolayer.radius);
  addQuantity("bulkEdge", &dropletPtr->bulk.radius);
  addQuantity("contactAngle", &dropletPtr->bulk.contactAngle);
  addQuantity("dropletHeight", &dropletPtr->bulk.height);
}

void FrameWriter::openFile(const char* _path) {
  joinPath(path, outputDir, _path);
  file = fopen(path, "w");
  if(options.verbose) {
    cout << "Opened file " << path << endl;
  }
}

void FrameWriter::closeFile() {
  fclose(file);
}

void FrameWriter::getDefaultFmt(char* final_fmt, double* dataPtr) {
  strcpy(final_fmt, "%15.6f");
}

void FrameWriter::getDefaultFmt(char* final_fmt, int* dataPtr) {
  strcpy(final_fmt, "%15d");
}

void FrameWriter::deleteFmtStrs() {
  vector<const char*>::iterator it;
  for(it=fmtArray.begin(); it<fmtArray.end(); it++) {
    delete *it;
  }
}

void FrameWriter::getQuantityStr(char* quantityStr, int i) {
  if(typeArray[i] == 'd') {
    sprintf(quantityStr, fmtArray[i], *((double*) dataPtrArray[i]));
  }
  else if(typeArray[i] == 'i') {
    sprintf(quantityStr, fmtArray[i], *((int*) dataPtrArray[i]));
  }
}

void FrameWriter::concatenateQuantityStrs() {
  char quantityStr[256];
  line[0] = 0;
  for(int i=0; i<numQuantities; i++) {
    getQuantityStr(quantityStr, i);
    strcat(line, quantityStr);
  }
}

// Write quantity names to first line of file
void FrameWriter::writeHeader() {
  char headerStr[256];
  line[0] = 0;
  for(int i=0; i<numQuantities; i++) {
    // Format quantity name (15 char wide)
    sprintf(headerStr, "%15s", quantityNameArray[i]);
    // Concatenate header strings
    strcat(line, headerStr);
  }
  // Write line to file
  fprintf(file, "%s\n", line);
}

void FrameWriter::writeFrame() {
  concatenateQuantityStrs();
  fprintf(file, "%s\n", line);
}
