
#include "Writers.h"

/////////////
// Writers //
/////////////

FrameWriter::FrameWriter() {
}

FrameWriter::setOutputDir(string path) {
  outputDir = path;
}

FrameWriter::openStream(string streamName, fmt="%15.6f") {
  stringstream pathStream;
  FILE* file;

  pathStream << outputDir << "/" << streamName << ".txt";
  file = fopen(pathStream.str(), "w");
  streams[streamName] = file;
}

FrameWriter::closeStream(string streamName) {
  close(streams[streamName]);
}

template <typename T>
FrameWriter::write(string streamName, T data) {
  fprintf(streams[streamName], fmts[streamName], data);
}

template <typename T>
FrameWriter::write(string streamName, vector<T> data, string fmt) {
  typename vector<T>::iterator it;
  for(it=data.begin(); it<data.end(); it++) {
    fprintf(streams[streamName], fmts[streamName], *it);
  }
}

