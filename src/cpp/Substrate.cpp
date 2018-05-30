#include "Substrate.h"

Substrate::Substrate(AtomArray &atomArray) {
  setContext(atomArray);
}

void Substrate::setContext(AtomArray &atomArray) {
  atomArrayPtr = &atomArray;
  simDataPtr = atomArrayPtr->simDataPtr;
  options = simDataPtr->options;
}
