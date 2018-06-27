#include "Writers.h"

/////////////
// Writers //
/////////////

WriterBase::WriterBase(SimData &simData) {
  numQuantities = 0;
  simDataPtr = &simData;
  options = simDataPtr->options;
  setOutputDir(options.outLoc.data());
  createDataDir();
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

// Create `outputDir`/`subdir`/`title`
void WriterBase::createDataDir() {
  char dataDir[256];
  // Create `outputDir`
  if(!dir_exists(outputDir)) {
    mkdir(outputDir, S_IRWXU);
  }

  // create `outputDir`/data
  joinPath(dataDir, outputDir, "data");
  if(!dir_exists(dataDir)) {
    mkdir(dataDir, S_IRWXU);
  }
}

void WriterBase::getFmtStr(char* fmt, double* dataPtr) {
  char lj[2];
  lj[0] = leftJustify ? '-' : '\0';
  lj[1] = 0;
  sprintf(fmt, "%%%s%d.%d%s", lj, COL_WIDTH, PRECISION, "f");
}

void WriterBase::getFmtStr(char* fmt, int* dataPtr) {
  char lj[2];
  lj[0] = leftJustify ? '-' : '\0';
  lj[1] = 0;
  sprintf(fmt, "%%%s%d%s", lj, COL_WIDTH, "d");
}
void WriterBase::getFmtStr(char* fmt, char* dataPtr) {
  char lj[2];
  lj[0] = leftJustify ? '-' : '\0';
  lj[1] = 0;
  sprintf(fmt, "%%%s%d%s", lj, COL_WIDTH, "s");
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

  leftJustify = true;

  setOutputQuantities();
  writeHeaders();

}

ScalarWriter::~ScalarWriter() {
  closeFiles();
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
  if(options.monolayer) {
    addQuantity("monoEdge", &dropletPtr->monolayer.radius);
  }
  if(options.fitCircle) {
    addQuantity("bulkEdge", &dropletPtr->bulk.radius);
    addQuantity("contactAngle", &dropletPtr->bulk.contactAngle);
    addQuantity("dropletHeight", &dropletPtr->bulk.height);
  }
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
  getFmtStr(fmtStr, headerStr);
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
    fprintf(files[i], "%s\n", quantityStr);
  }
}

//////////////////
// ArrayWriter //
//////////////////

ArrayWriter::ArrayWriter(DumpfileReader &dumpfileReader, Droplet &droplet) : WriterBase(*droplet.simDataPtr) {
  dumpfileReaderPtr = &dumpfileReader;
  dropletPtr = &droplet;

  leftJustify = true;

  setOutputQuantities();
}

ArrayWriter::~ArrayWriter() {
}

void ArrayWriter::setOutputQuantities() {
  addQuantity("boundaryPoints", dropletPtr->bulk.boundaryPointsArray, &dropletPtr->bulk.numPoints, dropletPtr->bulk.headers, 2);
}

// Get formatted string for quantity column
void ArrayWriter::getQuantityStr(char* quantityStr, int quantityNum, int i, int j) {
  double **doubleArr;
  int **intArr;

  // TODO: Replace 2 w/ numColumns
  // TODO: Make sure this is right
  if(typeArray[quantityNum] == 'd') {
    doubleArr = (double**) dataPtrArray[quantityNum];
    sprintf(quantityStr, fmtArray[quantityNum], doubleArr[j][i]);
  }
  else if(typeArray[quantityNum] == 'i') {
    intArr = (int**) dataPtrArray[quantityNum];
    sprintf(quantityStr, fmtArray[quantityNum], intArr[j][i]);
  }
}

void ArrayWriter::concatenateQuantityStrs(int quantityNum, int i) {
  char quantityStr[256];
  line[0] = 0;
  for(int j=0; j<numColumnsArray[quantityNum]; j++) {
    getQuantityStr(quantityStr, quantityNum, i, j);
    strcat(line, quantityStr);
  }
}

void ArrayWriter::createQuantityDir(const char* quantityName) {
  char quantityDir[256];
  char dataDir[256];

  // create `outputDir`/data
  joinPath(dataDir, outputDir, "data");
  if(!dir_exists(dataDir)) {
    mkdir(dataDir, S_IRWXU);
  }

  // create `outputDir`/data/`quantityName`
  joinPath(quantityDir, dataDir, quantityName);
  if(!dir_exists(quantityDir)) {
    mkdir(quantityDir, S_IRWXU);
  }
}

void ArrayWriter::getFilePath(char* filePath, const char* quantityName) {
  char filename[256];
  char dataDir[256];
  char quantityDir[256];

  joinPath(dataDir, outputDir, "data");
  joinPath(quantityDir, dataDir, quantityName);
  sprintf(filename, "%09d.txt", simDataPtr->framePtr->time);
  joinPath(filePath, quantityDir, filename);
}

// Write quantity names to first line of files
void ArrayWriter::writeHeader(FILE* file, int quantityNum) {
  char headerStr[256];
  char fmtStr[16];
  int numColumns = numColumnsArray[quantityNum];
  getFmtStr(fmtStr, headerStr);
  // Prepend # (comment line)
  for(int j=0; j<numColumns; j++) {
    sprintf(headerStr, "#%s", headersArray[quantityNum][j]);
    fprintf(file, fmtStr, headerStr);
  }
  fprintf(file, "\n");
}

void ArrayWriter::writeQuantityFile(int quantityNum) {
  char quantityStr[256];
  char filePath[265];
  int numRows = *lengthPtrArray[quantityNum];
  getFilePath(filePath, quantityNameArray[quantityNum]);
  FILE* file = fopen(filePath, "w");
  if(options.verbose) {
    cout << "Opened file '" << filePath << "' for array writing." << endl;
  }
  writeHeader(file, quantityNum);
  for(int i=0; i<numRows; i++) {
    concatenateQuantityStrs(quantityNum, i);
    fprintf(file, "%s\n", line);
  }
  fclose(file);
}

void ArrayWriter::writeFrame() {
  for(int quantityNum=0; quantityNum<numQuantities; quantityNum++) {
    writeQuantityFile(quantityNum);
  }
}
