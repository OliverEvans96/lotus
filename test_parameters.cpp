#include "Parameters.h"

using namespace std;

int main(int argc, char* argv[]) {
  CommandLineParser commandLineParser(argc, argv);
  commandLineParser.print();
  return 0;
}
