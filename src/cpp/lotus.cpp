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

  cout << "Initializing datafile" << endl;
  DatafileReader datafileReader(simData);
  cout << "Initializing atoms" << endl;
  AtomArray atoms(simData);
  cout << "Initializing dumpfile" << endl;
  DumpfileReader dumpfileReader(atoms);

  cout << "Initializing substrate" << endl;
  Substrate substrate(atoms, 5);
  cout << "Initializing droplet" << endl;
  Droplet droplet(atoms);

  cout << "Initializing figures" << endl;
  DensFigure densFigure("dens", droplet, substrate);
  DropletFigure dropletFigure("droplet", droplet);
  TanhFigure tanhFigure("tanh", droplet.monolayer.tanhFit);

  cout << "Initializing writers" << endl;
  ScalarWriter scalarWriter(dumpfileReader, droplet);
  ArrayWriter arrayWriter(dumpfileReader, droplet);

  cout << "Time loop" << endl;
  while(dumpfileReader.good()) {
    cout << "Reading frame" << endl;
    dumpfileReader.readFrame();

    if(options.substrate) {
      cout << "Filling substrate" << endl;
      substrate.fill(atoms);
    }

    cout << "Filling droplet" << endl;
    droplet.fill(atoms);

    if(options.monolayer) {
      cout << "Finding monolayer" << endl;
      droplet.findMonolayer();
      droplet.monolayer.fill(atoms);
      droplet.monolayer.calculateRadius();
    }

    cout << "Droplet calculations" << endl;
    droplet.dropletCalculations();

    cout << "Drawing figures" << endl;
    cout << "densFigure draw" << endl;
    densFigure.draw();
    cout << "densFigure save" << endl;
    densFigure.save();

    cout << "tanhFigure draw" << endl;
    tanhFigure.draw();
    cout << "tanhFigure save" << endl;
    tanhFigure.save();

    cout << "dropletFigure draw" << endl;
    dropletFigure.draw();
    cout << "dropletFigure save" << endl;
    dropletFigure.save();

    cout << "Writing data" << endl;
    scalarWriter.writeFrame();
    arrayWriter.writeFrame();
  }
}

