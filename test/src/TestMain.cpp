#include "Parameters.h"
#include "MDBase.h"
#include "Atoms.h"
#include "Readers.h"
#include "Droplet.h"
#include "Substrate.h"
#include "Visualization.h"
#include "catch.hpp"

using namespace std;

const int NUM_ATOMS = 513699;
const int NUM_STEPS = 8;
const int NUM_WATER = 1105;
const int STEPS_PER_FRAME = 3;
const double SUBSTRATE_MASS = 3130904.0;
const double DROPLET_MASS = 6563.56701333;

TEST_CASE("Readers", "[lotus]") {
  int argc = 2;
  char* argv[2] = {"test", "../test/data/test_config.yaml"};
  CommandLineParser commandLineParser(argc, argv);
  Options options = commandLineParser.options;
  REQUIRE(options.configPath == argv[1]);
  REQUIRE(options.liquidTypes[0] == 4);
  REQUIRE(options.liquidTypes[1] == 5);
  REQUIRE(options.solidTypes[0] == 1);
  REQUIRE(options.solidTypes[1] == 2);
  REQUIRE(options.solidTypes[2] == 3);
  REQUIRE(options.geometry == "spherical");
  REQUIRE(options.stepsPerFrame == STEPS_PER_FRAME);
  REQUIRE(options.dumpfile == "../test/data/20A_atom1_13-20");
  REQUIRE(options.datafile == "../test/data/lammps_noZperiod_3A.dat");
  REQUIRE(options.outLoc == "../test/run");
  REQUIRE(options.skipToEnd == false);
  REQUIRE(options.trackMonoAtoms == false);
  REQUIRE(options.saveImages == true);
  REQUIRE(options.plotHist == false);
  REQUIRE(options.plotDipole == false);
  REQUIRE(options.plotVr == false);
  REQUIRE(options.plotDensity == false);
  REQUIRE(options.plotAllTogether == false);
  // REQUIRE(options.verbose == false);
  REQUIRE(options.onlyFindInterface == false);
  SimData simData(options);
  DatafileReader datafileReader(simData);
  REQUIRE(simData.masses.size() == 5);
  REQUIRE(simData.masses[1] == 26.981540);
  REQUIRE(simData.masses[2] == 15.994915);
  REQUIRE(simData.masses[3] == 1.007825);
  REQUIRE(simData.masses[4] == 15.9994);
  REQUIRE(simData.masses[5] == 1.00794);
  REQUIRE(simData.numAtoms == NUM_ATOMS);
  REQUIRE(simData.waterBonds.size() == NUM_WATER);

  AtomArray atoms(simData);
  REQUIRE(atoms.numAtoms == NUM_ATOMS);

  DumpfileReader dumpfileReader(atoms);
  REQUIRE(simData.lastFrame.frameNum == 1);
  REQUIRE(simData.lastFrame.numSteps == 5);
  REQUIRE(simData.numSteps == NUM_STEPS);

  Substrate substrate(atoms, 5);
  Droplet droplet(atoms);

  char filename[100];

  cout << "*A" << endl;
  DensFigure densFigure("testDensFig", droplet, substrate);
  cout << "*B" << endl;

  DropletFigure dropletFigure("testHistFig", droplet);
  cout << "*C" << endl;

  TanhFigure tanhFigure("tanhFigure", droplet.monolayer.tanhFit);
  cout << "*D" << endl;

  // Time loop
  while(dumpfileReader.good()) {
    dumpfileReader.readFrame();
    substrate.fill(atoms);
    // REQUIRE(abs(substrate.getMass() - SUBSTRATE_MASS) < 1e-3*SUBSTRATE_MASS);
    droplet.fill(atoms);
    // REQUIRE(abs(droplet.getMass() - DROPLET_MASS) < 1e-3*DROPLET_MASS);
    // This is just mass in cylinder, not total
    // TODO: Fill monolayer
    cout << "----- FM -----" << endl;
    droplet.findMonolayer();
    cout << "----- FM-out -----" << endl;

    // TODO: Save image
    densFigure.draw();
    densFigure.save("dens.png");

    droplet.monolayer.fill(atoms);
    REQUIRE(droplet.monolayer.hMono->GetEntries() > 0);
    REQUIRE(droplet.monolayer.zlim[1] - droplet.monolayer.zlim[0] > 0);
    cout << "DATA?" << endl;
    TH1D *h = droplet.monolayer.tanhFit.hTanh;
    int n = h->GetNbinsX();
    double *x = new double[n];
    double *y = new double[n];
    for(int i=1; i<h->GetNbinsX(); i++) {
      cout << h->GetBinCenter(i) << " - " << h->GetBinContent(i) << endl;
      x[i-1] = h->GetBinCenter(i);
      y[i-1] = h->GetBinContent(i);
    }

    cout << "TEST FIT" << endl;
    TGraph *g = new TGraph(n, x, y);

    cout << "~A" << endl;
    TH1D *h1 = new TH1D("test", "test", 10, 0, 50);
    cout << "~B" << endl;
    h1->Fill(1.0);
    h1->Fill(12.0);
    h1->Fill(13.0);
    h1->Fill(29.0);
    h1->Fill(36.0);
    h1->Fill(16.0);
    h1->Fill(27.0);
    h1->Fill(54.0);
    h1->Fill(18.0);
    h1->Fill(32.0);
    cout << "f1" << endl;

    h1->Fit("gaus");
    cout << "gh" << endl;

    delete h1;

    // cout << "-=A" << endl;
    // g->Fit("gaus");
    cout << "-=B" << endl;
    g->Fit(droplet.monolayer.tanhFit.fTanh, "W");
    // droplet.monolayer.tanhFit.hTanh->Fit(droplet.monolayer.tanhFit.fTanh, "W");

    cout << "deallocate" << endl;
    delete g;

    delete [] x;
    delete [] y;
    cout << "END TEST" << endl;

    droplet.monolayer.calculateRadius();

    tanhFigure.draw();
    tanhFigure.save("mono_tanh.png");

    REQUIRE(droplet.monolayer.radius > 0);

    droplet.dropletCalculations();

    // TODO: Draw droplet figure
    dropletFigure.draw();
    dropletFigure.save("droplet.png");


    REQUIRE(droplet.bulk.gCirclePoints->GetN() > 0);
    REQUIRE(droplet.bulk.circle.GetNumPoints() > 0);

    REQUIRE(droplet.bulk.circle.intersected);
    REQUIRE(droplet.bulk.height > 0);
    REQUIRE(droplet.bulk.radius > 0);
    REQUIRE(droplet.bulk.contactAngle >= -1);
    REQUIRE(droplet.bulk.contactAngle <= 1);
    // TODO: Write frame quantities to file
    //       - rm, rb, ca, h, circle points
    // TODO: Add option to enable/disable circle fit

    // TODO: Clean up print statements

    // REQUIRE(file_exists("out.png"));

    // sprintf(filename, "substrate_density%d.png", dumpfileReader.frameNum);
    // substrate.plotDensity(filename);

    // sprintf(filename, "droplet_density%d.png", dumpfileReader.frameNum);
    // droplet.plotDensity(filename);

  }
}

