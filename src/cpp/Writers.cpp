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
  sprintf(fmt, "%%%s%d.%d%s", lj, options.outputColWidth, options.outputPrecision, "f");
}

void WriterBase::getFmtStr(char* fmt, int* dataPtr) {
  char lj[2];
  lj[0] = leftJustify ? '-' : '\0';
  lj[1] = 0;
  sprintf(fmt, "%%%s%d%s", lj, options.outputColWidth, "d");
}
void WriterBase::getFmtStr(char* fmt, const char* dataPtr) {
  char lj[2];
  lj[0] = leftJustify ? '-' : '\0';
  lj[1] = 0;
  sprintf(fmt, "%%%s%d%s", lj, options.outputColWidth, "s");
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

/**
   Save appropriate pointers,
   left justify scalar data,
   set output quantities,
   and write appropriate headers.
*/
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

/// Only close files which are actually open.
void ScalarWriter::closeFiles() {
  for(int i=0; i<numQuantities; i++) {
    if(files[i] != NULL) {
      fclose(files[i]);
    }
  }
}

/**
   A quantity name and pointer to the quantity are given for each.

   @see #addQuantity
*/
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

/**
   Since #dataPtrArray is a vector of void pointers,
   this function checks #typeArray in order to know
   how to appropriately dereference the pointer.

   Also, values larger than @p 1e20 are replaced with "INF"
   (formatted as the quantity would have been).

   @param[out] quantityStr    Final formatted quantity string.
   @param[in]  i              Index of the quantity.
*/
void ScalarWriter::getQuantityStr(char* quantityStr, int i) {
  double doubleVal;
  int intVal;

  // Don't print quantities larger than this
  double maxVal = 1e20;
  char fmt[16];
  const char infStr[16] = "INF";

  if(typeArray[i] == 'd') {
    doubleVal = *((double*) dataPtrArray[i]);
    if(doubleVal < maxVal) {
      sprintf(quantityStr, fmtArray[i], doubleVal);
    }
    else {
      getFmtStr(fmt, infStr);
      sprintf(quantityStr, fmt, infStr);
    }
  }
  else if(typeArray[i] == 'i') {
    intVal = *((int*) dataPtrArray[i]);
    if(intVal < maxVal) {
      sprintf(quantityStr, fmtArray[i], intVal);
    }
    else {
      getFmtStr(fmt, infStr);
      sprintf(quantityStr, fmt, infStr);
    }
  }
}

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
    fflush(files[i]);
  }
}

/// `FILE*`s are flushed after each write.
void ScalarWriter::writeFrame() {
  char quantityStr[256];
  for(int i=0; i<numQuantities; i++) {
    getQuantityStr(quantityStr, i);
    fprintf(files[i], "%s\n", quantityStr);
    fflush(files[i]);
  }
}

//////////////////
// ArrayWriter //
//////////////////

/**
   Save appropriate pointers,
   left justify scalar data,
   and set output quantities.
*/
ArrayWriter::ArrayWriter(DumpfileReader &dumpfileReader, Droplet &droplet) : WriterBase(*droplet.simDataPtr) {
  dumpfileReaderPtr = &dumpfileReader;
  dropletPtr = &droplet;

  leftJustify = true;

  setOutputQuantities();
}

ArrayWriter::~ArrayWriter() {
}

/**
   Each quantity requires:
   - A quantity name, a pointer to the data array,
   - a pointer to the number of rows,
   - a pointer to the array of header strings,
   - the number of columns.

   Currently, only @p boundaryPoints (CircularBulk::boundaryPointsArray)
   is written in this manner.

   @see #addQuantity
*/
void ArrayWriter::setOutputQuantities() {
  addQuantity(
      "boundaryPoints",
      dropletPtr->bulk.boundaryPointsArray,
      &dropletPtr->bulk.numPoints,
      dropletPtr->bulk.headers,
      2
  );
}

/**
   Since #dataPtrArray is a vector of void pointers,
   this function checks #typeArray in order to know
   how to appropriately dereference the pointer.

   Also, values larger than @p 1e20 are replaced with "INF"
   (formatted as the quantity would have been).

   @param[out]  quantityStr    Final formatted string for this element.
   @param[in]   quantityNum    Index of the quantity.
   @param[in]   i              Row number of this element.
   @param[in]   j              Column number of this element.
*/
void ArrayWriter::getQuantityStr(char* quantityStr, int quantityNum, int i, int j) {
  double **doubleArr;
  int **intArr;

  double doubleVal;
  int intVal;

  // Don't print quantities larger than this
  double maxVal = 1e20;
  char fmt[16];
  const char infStr[16] = "INF";

  if(typeArray[quantityNum] == 'd') {
    doubleArr = (double**) dataPtrArray[quantityNum];
    doubleVal = doubleArr[j][i];
    if(doubleVal < maxVal) {
      sprintf(quantityStr, fmtArray[quantityNum], doubleArr[j][i]);
    }
    else {
      getFmtStr(fmt, infStr);
      sprintf(quantityStr, fmt, infStr);
    }
  }
  else if(typeArray[quantityNum] == 'i') {
    intArr = (int**) dataPtrArray[quantityNum];
    intVal = intArr[j][i];
    if(intVal < maxVal) {
      sprintf(quantityStr, fmtArray[quantityNum], intArr[j][i]);
    }
    else {
      getFmtStr(fmt, infStr);
      sprintf(quantityStr, fmt, infStr);
    }
  }
}

/**
   Final formatted quantity string is saved to #quantityStr.

   @param quantityNum   Index of the quantity.
   @param i             Row number of the quantity
*/
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

/**
   @param[out] filePath       Final path to quantity file for this timestep.
   @param[in]  quantityName   Name of this quantity
*/
void ArrayWriter::getFilePath(char* filePath, const char* quantityName) {
  char filename[256];
  char dataDir[256];
  char quantityDir[256];

  joinPath(dataDir, outputDir, "data");
  joinPath(quantityDir, dataDir, quantityName);
  sprintf(filename, "%09d.txt", simDataPtr->framePtr->time);
  joinPath(filePath, quantityDir, filename);
}

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
  fflush(file);
  fclose(file);
}

void ArrayWriter::writeFrame() {
  for(int quantityNum=0; quantityNum<numQuantities; quantityNum++) {
    writeQuantityFile(quantityNum);
  }
}
