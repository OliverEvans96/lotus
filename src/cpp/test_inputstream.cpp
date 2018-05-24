#include "Utils.h"
#include <iostream>

using namespace std;

int main() {
  cout << "1" << endl;
  InputStream goodStream("centroid.txt");
  cout << "2" << endl;
  InputStream badStream("fakefile.txt");
  cout << "3" << endl;
  return 0;
}
