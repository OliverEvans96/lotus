#include "Parameters.h"

///////////////////////////
// High level parameters //
///////////////////////////

//Options::Options(char* optionsFile) {
//}

// Read a YAML file into a map of string vectors.
// Does not support nested YAML mappings.
StrVecMap Options::parseYaml(const char* filename) {
    FILE *config_file;
    yaml_parser_t parser;
    yaml_document_t document;
    yaml_node_t *node, *node1, *key, *value;
    yaml_node_item_t *item;
    yaml_node_pair_t *pair;
    string str;
    vector<string> vec;
    StrVecMap yamlMap;

    config_file = fopen(filename, "rb");

    yaml_parser_initialize(&parser);
    yaml_parser_set_input_file(&parser, config_file);
    yaml_parser_load(&parser, &document);

    node = yaml_document_get_root_node(&document);
    assert(node->type == YAMLMAPPING_NODE);

    // Mapping
    for(pair=node->data.mapping.pairs.start; pair<node->data.mapping.pairs.top; pair++) {
        key = yaml_document_get_node(&document, pair->key);
        value = yaml_document_get_node(&document, pair->value);
        // Scalar
        if(value->type == YAML_SCALAR_NODE) {
            vec.push_back((char*)value->data.scalar.value);
        }
        // Sequence
        else if(value->type == YAML_SEQUENCE_NODE) {
            for(item=value->data.sequence.items.start;
                    item<value->data.sequence.items.top; item++) {
                node1 = yaml_document_get_node(&document, *item);
                vec.push_back((char*)node1->data.scalar.value);
            }
        }
        yamlMap[(char*)key->data.scalar.value] = vec;
        vec.clear();
    }

    yaml_document_delete(&document);
    yaml_parser_delete(&parser);
    fclose(config_file);

    return yamlMap;
}

// Determine whether key is present in map
bool Options::mapHasKey(StrVecMap yamlMap, string key) {
  map<string, vector<string> >::iterator it;
  it = yamlMap.find(key);
  return (it != yamlMap.end())
}

void Options::printYamlMap(StrVecMap yamlMap) {
    string key, val;
    vector<string> vec;
    StrVecMap::iterator map_it;
    vector<string>::iterator vec_it;

    for(map_it=yamlMap.begin(); map_it != yamlMap.end(); map_it++) {
        key = map_it -> first;
        vec = map_it -> second;
        printf("%s: [", key.data());
        for(vec_it=vec.begin(); vec_it != vec.end(); vec_it++) {
           val = *vec_it;
           printf("\"%s\", ", val.data());
        }
        printf("]\n");
    }
}

// Convert strings to appropriate types
void Options::fromString(string optionString, bool &option) {
  stringstream ss;
  ss.str(optionString);
  ss >> boolalpha >> option;
}
void Options::fromString(string optionString, int &option) {
  option = atoi(optionString.data());
}
void Options::fromString(string optionString, double &option) {
  option = atof(optionString.data());
}
void Options::fromString(string optionString, string &option) {
  option = optionString.data();
}

// Parse scalar option, assuming key is present
template <typename T>
void Options::unsafeParseOption(string optionName, T &option) {
  string optionString;
  optionString = yamlMap[optionName][0];
  fromString(optionName, option);
}

// Parse vector option, assuming key is present
template <typename T>
void Options::unsafeParseOption(string optionName, vector<T> &optionVec) {
  vector<string> strVec;
  vector<string>::iterator it;
  T option;

  strVec = yamlMap[optionName];
  for(it=strVec.begin(); it<strVec.end(); it++) {
    fromString(*it, option);
    optionVec.push_back(option);
  }
}

// Safely parse option, falling back on default value
template <typename T>
void Options::parseDefaultOption(string optionName, T &option, T &defaultValue) {
  if(mapHasKey(yamlMap, optionName)) {
    unsafeParseOption(optionName, option);
  }
  else {
    option = defaultValue;
  }
}
// Safely parse option, throwing exception if unspecified.
template <typename T>
void Options::parseRequiredOption(string optionName, T &option) {
  stringstream err;

  if(mapHasKey(yamlMap, optionName)) {
    unsafeParseOption(optionName, option);
  }
  else {
    err << "'" << optionName << "' is required in YAML config (" << optionsFile << ").";
    throw ss.str().data();
  }
}

void Options::readConfig(string configFile) {
  yamlMap = parseYaml(configFile);

  parseRequiredOption("liquidTypes", liquidTypes);
  parseRequiredOption("solidTypes", solidTypes);
  parseRequiredOption("inLoc", inLoc);
  parseRequiredOption("outLoc", outLoc);

  parseDefaultOption("skipToEnd", skipToEnd, false);
  parseDefaultOption("trackMonoAtoms", trackMonoAtoms, false);
  parseDefaultOption("saveImages", saveImages, false);
  parseDefaultOption("plotHist", plotHist, false);
  parseDefaultOption("plotDipole", plotDipole, false);
  parseDefaultOption("plotVr", plotVr, false);
  parseDefaultOption("plotDensity", plotDensity, false);
  parseDefaultOption("plotAllTogether", plotAllTogether, false);
  parseDefaultOption("debugOutput", debugOutput, false);
  parseDefaultOption("onlyFindInterface", onlyFindInterface, false);
}

void Options::printOptions() {
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
  optionsFile = argv[1];
}

void CommandLineParser::print() {
  cout << "Command Line Arguments: " << endl;
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
