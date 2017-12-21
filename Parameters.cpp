#include "Parameters.h"

///////////////////////////
// High level parameters //
///////////////////////////

//Options::Options(char* optionsFile) {
//}

void Options::readOptions(string optionsFile) {
  // Format:
  // optionName bool

  ifstream optionsStream(optionsFile);
  string firstLine, option;
  bool flag;

  // Ignore first line
  getline(optionsStream, firstLine);

  while(!optionsStream.eof()) {
    optionsStream >> option >> flag;
    optionsStream.ignore(256, '\n');

    if(option=="skipToEnd") {
      skipToEnd = flag;
    }
    else if(option=="trackMonoAtoms") {
      trackMonoAtoms = flag;
    }
    else if(option=="saveImages") {
      saveImages = flag;
    }
    else if(option=="plotHist") {
      plotHist = flag;
    }
    else if(option=="plotDipole") {
      plotDipole = flag;
    }
    else if(option=="plotVr") {
      plotVr = flag;
    }
    else if(option=="plotDensity") {
      plotDensity = flag;
    }
    else if(option=="plotAllTogether") {
      plotAllTogether = flag;
    }
    else if(option=="debugOutput") {
      debugOutput = flag;
    }
    else if(option=="onlyFindInterface") {
      onlyFindInterface = flag;
    }
  }
}

void Options::print() {
  cout << "Options:" << endl;
  cout << "skipToEnd = " << skipToEnd << endl;
  cout << "trackMonoAtoms = " << trackMonoAtoms << endl;
  cout << "saveImages = " << saveImages << endl;
  cout << "plotHist = " << plotHist << endl;
  cout << "plotDipole = " << plotDipole << endl;
  cout << "plotVr = " << plotVr << endl;
  cout << "plotDensity = " << plotDensity << endl;
  cout << "plotAllTogether = " << plotAllTogether << endl;
  cout << "debugOutput = " << debugOutput << endl;
  cout << "onlyFindInterface = " << onlyFindInterface << endl;
}

CommandLineParser::CommandLineParser(int argc, char* argv[]) {
  //Command line arguments: inLoc,outLoc

  /*
    cout << "argv: " << endl;
    for(int i=0;i<argc;i++)
    cout << argv[i] << endl;
    cout << endl;
  */

  cout << "Received " << argc-1 << " command line arguments." << endl;
  parseArgs(argc, argv);
  options.readOptions(optionsFile);
}

void CommandLineParser::parseArgs(int argc, char* argv[]) {
  inLoc = argv[1];
  outLoc = argv[2];
  optionsFile = argv[3];
}

void CommandLineParser::print() {
  cout << "Command Line Arguments: " << endl;
  cout << "inLoc = " << inLoc << endl;
  cout << "outLoc = " << outLoc << endl;
  cout << "optionsFile = " << optionsFile << endl;
  options.print();
}


// //Whether to run whole analysis or only find monolayer
// bool onlyFindInterface;
// cout << "argc=" << argc << endl;
// if(argc<4)
//   onlyFindInterface=false;
//  else
//    onlyFindInterface=(strcmp(argv[3],"mono")==0);


//While reading files, we don't yet know where the monolayer will be defined for this frame, so we save an atom's coordinates if it's in a broad range, then eliminate extras later.
//double potentialMonoLimits[2]={14.5,15.5};
//int nMonoAtoms=0;
//int nPotentialMonoAtoms=0;
//vector<double> potentialMonoR;
//vector<double> potentialMonoZ;

//Open monoLimits file
/*
  fstream interfaceFile;
  if(!onlyFindInterface)
  {
  interfaceFile.open("../mono_limits.txt",ios::in);
  }
  else
  interfaceFile.open("../mono_limits.txt",ios::out);


  if (interfaceFile.good())
  cout << "Successfully opened " << "mono_limits.txt" << endl;
  else
  cout << "Failed to open " << "mono_limits.txt" << endl;

  if(!onlyFindInterface)
  {
  interfaceFile >> monoLimits[0] >> monoLimits[1];
  cout << "monoLimits=" << monoLimits[0] << " " << monoLimits[1] << endl;
  }
*/



//AnalysisParameters::AnalysisParameters()
