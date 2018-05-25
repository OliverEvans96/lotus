
#include "Writers.h"

/////////////
// Writers //
/////////////

void FrameWriter::setOutputDir(string path) {
  outputDir = path;
}

void FrameWriter::openStream(string streamName, string fmt="%15.6f") {
  stringstream pathStream;

  pathStream << outputDir << "/" << streamName << ".txt";
  streams[streamName]  = fopen(pathStream.str().data(), "w");
  fmts[streamName] = fmt;
}

void FrameWriter::closeStream(string streamName) {
  fclose(streams[streamName]);
}

template <typename T>
void FrameWriter::write(string streamName, T data) {
  fprintf(streams[streamName], fmts[streamName], data);
}

template <typename T>
void FrameWriter::write(string streamName, vector<T> data, string fmt) {
  typename vector<T>::iterator it;
  for(it=data.begin(); it<data.end(); it++) {
    fprintf(streams[streamName], fmts[streamName], *it);
  }
}

