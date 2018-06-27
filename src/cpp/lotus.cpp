#include "Parameters.h"
#include "MDBase.h"
#include "Atoms.h"
#include "Readers.h"
#include "Writers.h"
#include "Droplet.h"
#include "Substrate.h"
#include "Visualization.h"

using namespace std;

int main(int argc, const char** argv) {
  CommandLineParser commandLineParser(argc, argv);
  Options options = commandLineParser.options;
  SimData simData(options);

  DatafileReader datafileReader(simData);
  AtomArray atoms(simData);
  DumpfileReader dumpfileReader(atoms);

  Substrate substrate(atoms, 5);
  Droplet droplet(atoms);

  DensFigure densFigure("dens", droplet, substrate);
  DropletFigure dropletFigure("droplet", droplet);
  TanhFigure tanhFigure("tanh", droplet.monolayer.tanhFit);

  ScalarWriter scalarWriter(dumpfileReader, droplet);
  ArrayWriter arrayWriter(dumpfileReader, droplet);

  // Time loop
  while(dumpfileReader.good()) {
    dumpfileReader.readFrame();

    if(options.substrate) {
      substrate.fill(atoms);
    }

    droplet.fill(atoms);

    if(options.monolayer) {
      droplet.findMonolayer();
      droplet.monolayer.fill(atoms);
      droplet.monolayer.calculateRadius();
    }

    droplet.dropletCalculations();

    densFigure.draw();
    densFigure.save();

    tanhFigure.draw();
    tanhFigure.save();

    dropletFigure.draw();
    dropletFigure.save();

    scalarWriter.writeFrame();
    arrayWriter.writeFrame();
  }
}

