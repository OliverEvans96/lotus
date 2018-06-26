
#include "Writers.h"

/////////////
// Writers //
/////////////

FrameWriter::FrameWriter(SimData &simData) {
  simDataPtr = &simData;
  options = simDataPtr->options;
  setOutputDir(options.outLoc.data());
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

void FrameWriter::stripTrailingSlash(char* strippedPath, char* path) {
  int len = strlen(path);

  strcpy(strippedPath, path);
  if(strippedPath[len-1] == '/') {
    // Replace with null terminator
    strippedPath[len-1] = 0;
  }
}

void FrameWriter::joinPath(char* path, char* prefix, char* suffix) {
  char strippedPath[256];
  // Remove trailing slash from prefix if present
  stripTrailingSlash(strippedPath, prefix);
  sprintf(path, "%s/%s", strippedPath, suffix);
}

void FrameWriter::setOutputDir(char* path) {
  strcpy(outputDir, path);
  if(options.verbose) {
    cout << "Set output dir '" << outputDir << "'" << endl;
  }
}

void FrameWriter::openFile(char* _path) {
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
  strcpy(final_fmt, "%10d");
}

void FrameWriter::deleteFmtStrs() {
  vector<const char*>::iterator it;
  for(it=fmtArray.begin(); it<fmtArray.end(); it++) {
    delete *it;
  }
}

void FrameWriter::getQuantityStr(char* quantityStr, int i) {
  cout << "t[" << i << "] = '" << typeArray[i] << "'" << endl;
  cout << "f[" << i << "] = '" << fmtArray[i] << "'" << endl;
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
    cout << "q=" << quantityNameArray[i] << endl;
    getQuantityStr(quantityStr, i);
    cout << "qs='" << quantityStr << "'" << endl;
    strcat(line, quantityStr);
  }
}

// TODO: Align columns or write to separate files
void FrameWriter::writeHeader() {
  char quantityStr[256];
  line[0] = 0;
  for(int i=0; i<numQuantities; i++) {
    strcat(line, quantityNameArray[i]);
    if(i<numQuantities-1)
      strcat(line, " ");
  }

  fprintf(file, "%s\n", line);
}

void FrameWriter::writeFrame() {
  concatenateQuantityStrs();
  fprintf(file, "%s\n", line);
  printf("Line = '%s'\n", line);
}
