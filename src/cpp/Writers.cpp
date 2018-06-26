
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
  strcpy(final_fmt, "%-15.6f");
}

void WriterBase::getDefaultFmt(char* final_fmt, int* dataPtr) {
  strcpy(final_fmt, "%-15d");
}
void WriterBase::getDefaultFmt(char* final_fmt, char* dataPtr) {
  strcpy(final_fmt, "%-15s");
}


void WriterBase::deleteFmtStrs() {
  vector<const char*>::iterator it;
  for(it=fmtArray.begin(); it<fmtArray.end(); it++) {
    delete *it;
  }
}


//////////////////
// ScalarWriter //
//////////////////

ScalarWriter::ScalarWriter(DumpfileReader &dumpfileReader, Droplet &droplet) : WriterBase(*droplet.simDataPtr) {
  dumpfileReaderPtr = &dumpfileReader;
  dropletPtr = &droplet;

  createDataDir();
  setOutputQuantities();
  writeHeaders();
}

ScalarWriter::~ScalarWriter() {
  closeFiles();
}

// Create `outputDir`/`subdir`/`title`
void ScalarWriter::createDataDir() {
  char dataDir[256];
  joinPath(dataDir, outputDir, "data");
  if(!dir_exists(dataDir)) {
    mkdir(dataDir, S_IRWXU);
  }
}

FILE* ScalarWriter::openFileBase(const char* filename) {
  char dataDir[256];
  joinPath(dataDir, outputDir, "data");
  joinPath(path, dataDir, filename);
  FILE* file = fopen(path, "w");
  if(options.verbose) {
    cout << "Opened file '" << path << "' for scalar writing." << endl;
  }
  return file;
}

void ScalarWriter::openFile(const char* quantityName) {
  char filename[256];
  sprintf(filename, "%s.txt", quantityName);
  files.push_back(openFileBase(filename));
}

void ScalarWriter::closeFiles() {
  for(int i=0; i<numQuantities; i++) {
    // Only close if file is actually open
    if(files[i] != NULL) {
      fclose(files[i]);
    }
  }
}

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

// Write quantity names to first line of files
void ScalarWriter::writeHeaders() {
  char headerStr[256];
  char fmtStr[16];
  getDefaultFmt(fmtStr, headerStr);
  strcat(fmtStr, "\n");
  for(int i=0; i<numQuantities; i++) {
    // Prepend # (comment line)
    sprintf(headerStr, "#%s", quantityNameArray[i]);
    // Format quantity name
    fprintf(files[i], fmtStr, headerStr);
  }
}

void ScalarWriter::writeFrame() {
  char quantityStr[256];
  for(int i=0; i<numQuantities; i++) {
    getQuantityStr(quantityStr, i);
    // Format quantity name (15 char wide)
    fprintf(files[i], "%15s\n", quantityStr);
  }
}

//////////////////
// VectorWriter //
//////////////////

void VectorWriter::getQuantityStr(char* quantityStr, int i) {
}

void VectorWriter::concatenateQuantityStrs() {
  char quantityStr[256];
  line[0] = 0;
  for(int i=0; i<numQuantities; i++) {
    getQuantityStr(quantityStr, i);
    strcat(line, quantityStr);
  }
}
