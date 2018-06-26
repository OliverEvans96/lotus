
#include "Writers.h"

/////////////
// Writers //
/////////////

WriterBase::WriterBase(SimData &simData) {
  simDataPtr = &simData;
  options = simDataPtr->options;
  setOutputDir(options.outLoc.data());
}

WriterBase::~WriterBase() {
  closeFile();
  deleteFmtStrs();
}

void WriterBase::storeType(double *x) {
  typeArray.push_back('d');
}
void WriterBase::storeType(int *x) {
  typeArray.push_back('i');
}

void WriterBase::setOutputDir(const char* path) {
  strcpy(outputDir, path);
  if(options.verbose) {
    cout << "Set output dir '" << outputDir << "'" << endl;
  }
}

void WriterBase::getDefaultFmt(char* final_fmt, double* dataPtr) {
  strcpy(final_fmt, "%15.6f");
}

void WriterBase::getDefaultFmt(char* final_fmt, int* dataPtr) {
  strcpy(final_fmt, "%15d");
}

void WriterBase::deleteFmtStrs() {
  vector<const char*>::iterator it;
  for(it=fmtArray.begin(); it<fmtArray.end(); it++) {
    delete *it;
  }
}

void WriterBase::openFile(const char* _path) {
  joinPath(path, outputDir, _path);
  file = fopen(path, "w");
  if(options.verbose) {
    cout << "Opened file " << path << endl;
  }
}

void WriterBase::closeFile() {
  fclose(file);
}


//////////////////
// ScalarWriter //
//////////////////

ScalarWriter::ScalarWriter(const char* filename, DumpfileReader &dumpfileReader, Droplet &droplet) : WriterBase(*droplet.simDataPtr) {
  dumpfileReaderPtr = &dumpfileReader;
  dropletPtr = &droplet;

  openFile(filename);

  setOutputQuantities();
  writeHeader();
}

ScalarWriter::~ScalarWriter() {}

// Choose which quantities are written to the main results file
// A quantity name and pointer to the quantity are given for each
void ScalarWriter::setOutputQuantities() {
  addQuantity("timestep", &dumpfileReaderPtr->frameReader.frame.time);
  addQuantity("monoEdge", &dropletPtr->monolayer.radius);
  addQuantity("bulkEdge", &dropletPtr->bulk.radius);
  addQuantity("contactAngle", &dropletPtr->bulk.contactAngle);
  addQuantity("dropletHeight", &dropletPtr->bulk.height);
}

void ScalarWriter::getQuantityStr(char* quantityStr, int i) {
  if(typeArray[i] == 'd') {
    sprintf(quantityStr, fmtArray[i], *((double*) dataPtrArray[i]));
  }
  else if(typeArray[i] == 'i') {
    sprintf(quantityStr, fmtArray[i], *((int*) dataPtrArray[i]));
  }
}

void ScalarWriter::concatenateQuantityStrs() {
  char quantityStr[256];
  line[0] = 0;
  for(int i=0; i<numQuantities; i++) {
    getQuantityStr(quantityStr, i);
    strcat(line, quantityStr);
  }
}

// Write quantity names to first line of file
void ScalarWriter::writeHeader() {
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

void ScalarWriter::writeFrame() {
  concatenateQuantityStrs();
  fprintf(file, "%s\n", line);
}
