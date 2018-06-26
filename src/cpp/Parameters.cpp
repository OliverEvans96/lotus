#include "Parameters.h"

///////////////////////////
// High level parameters //
///////////////////////////

//Options::Options(char* configPath) {
//}

// Read a YAML file into a map of string vectors.
// Does not support nested YAML mappings.
StrVecMap Options::parseYaml(const char* configPath) {
    FILE *configFile;
    yaml_parser_t parser;
    yaml_document_t document;
    yaml_node_t *node, *node1, *key, *value;
    yaml_node_item_t *item;
    yaml_node_pair_t *pair;
    string str;
    vector<string> vec;
    StrVecMap yamlMap;

    configFile = fopen(configPath, "rb");

    yaml_parser_initialize(&parser);
    yaml_parser_set_input_file(&parser, configFile);
    yaml_parser_load(&parser, &document);

    node = yaml_document_get_root_node(&document);
    assert(node->type == YAML_MAPPING_NODE);

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
    fclose(configFile);

    return yamlMap;
}

// Determine whether key is present in map
bool Options::mapHasKey(StrVecMap yamlMap, string key) {
  map<string, vector<string> >::iterator it;
  it = yamlMap.find(key);
  return (it != yamlMap.end());
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
  fromString(optionString, option);
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
void Options::parseDefaultOption(string optionName, T &option, T defaultValue) {
  if(mapHasKey(yamlMap, optionName)) {
    unsafeParseOption(optionName, option);
  }
  else {
    option = defaultValue;
  }
}
// Safely parse option, terminating if unspecified.
template <typename T>
void Options::parseRequiredOption(string optionName, T &option) {
  if(mapHasKey(yamlMap, optionName)) {
    unsafeParseOption(optionName, option);
  }
  else {
    cout << "'" << optionName << "' is required in YAML config (" << configPath << ")." << endl;
    exit(1);
  }
}

void Options::readConfig(const char* _configPath) {
  strcpy(configPath, _configPath);
  yamlMap = parseYaml(configPath);

  parseRequiredOption("liquidTypes", liquidTypes);
  parseRequiredOption("solidTypes", solidTypes);
  parseRequiredOption("dumpfile", dumpfile);
  parseRequiredOption("datafile", datafile);
  parseRequiredOption("outLoc", outLoc);

  parseDefaultOption("geometry", geometry, (string)"spherical");
  parseDefaultOption("stepsPerFrame", stepsPerFrame, 5);
  parseDefaultOption("waterBondType", waterBondType, 2);
  parseDefaultOption("skipToEnd", skipToEnd, false);
  parseDefaultOption("trackMonoAtoms", trackMonoAtoms, false);
  parseDefaultOption("saveImages", saveImages, false);
  parseDefaultOption("saveROOT", saveROOT, false);
  parseDefaultOption("plotHist", plotHist, false);
  parseDefaultOption("plotDipole", plotDipole, false);
  parseDefaultOption("plotVr", plotVr, false);
  parseDefaultOption("plotDensity", plotDensity, false);
  parseDefaultOption("plotAllTogether", plotAllTogether, false);
  parseDefaultOption("verbose", verbose, false);
  parseDefaultOption("onlyFindInterface", onlyFindInterface, false);
  parseDefaultOption("dz", dz, 1.0);
  parseDefaultOption("dv", dv, 250.0);
  parseDefaultOption("dens_min", dens_min, 0.0);
  parseDefaultOption("dens_max", dens_max, 5.0);
  parseDefaultOption("plot_rmax", plot_rmax, 150.0);
  parseDefaultOption("plot_zmax", plot_zmax, 100.0);
  parseDefaultOption("plot_aspect", plot_aspect, 1.0);
  parseDefaultOption("plot_width", plot_width, 800);
  parseDefaultOption("expectedLiquidDensity", expectedLiquidDensity, 1.0);

  // Set ROOT verbosity
  if(verbose) {
    gErrorIgnoreLevel = kInfo;
  }
  else {
    gErrorIgnoreLevel = kWarning;
  }
}


template <typename T>
void Options::printOption(string optionName, T option) {
  cout << optionName << ": " << option << endl;
}

template <typename T>
void Options::printOption(string optionName, vector<T> optionVec) {
  typename vector<T>::iterator it;
  it = optionVec.begin();
  cout << optionName << ": [";
  for(int i=0; i<optionVec.size()-1; i++) {
    cout << *it++ << ", ";
  }
  cout << *it << "]" << endl;
}

void Options::printOptions() {
  printOption("liquidTypes", liquidTypes);
  printOption("solidTypes", solidTypes);
  printOption("dumpfile", dumpfile);
  printOption("datafile", datafile);
  printOption("outLoc", outLoc);
  printOption("waterBondType", waterBondType);
  printOption("geometry", geometry);
  printOption("stepsPerFrame", stepsPerFrame);
  printOption("skipToEnd", skipToEnd);
  printOption("trackMonoAtoms", trackMonoAtoms);
  printOption("saveImages", saveImages);
  printOption("saveROOT", saveROOT);
  printOption("plotHist", plotHist);
  printOption("plotDipole", plotDipole);
  printOption("plotVr", plotVr);
  printOption("plotDensity", plotDensity);
  printOption("plotAllTogether", plotAllTogether);
  printOption("verbose", verbose);
  printOption("onlyFindInterface", onlyFindInterface);
  printOption("dens_min", dens_min);
  printOption("dens_max", dens_max);
  printOption("plot_rmax", plot_rmax);
  printOption("plot_zmax", plot_zmax);
  printOption("plot_aspect", plot_aspect);
  printOption("dz", dz);
  printOption("dv", dv);
  printOption("plot_rmax", plot_rmax);
  printOption("plot_zmax", plot_zmax);
  printOption("plot_aspect", plot_aspect);
  printOption("plot_width", plot_width);
  printOption("expectedLiquidDensity", expectedLiquidDensity);
}

CommandLineParser::CommandLineParser(int argc, const char* argv[]) {
  parseArgs(argc, argv);
  options.readConfig(configPath);
  if(options.verbose) {
    print();
  }
}

void CommandLineParser::parseArgs(int argc, const char* argv[]) {
  if(argc == 2) {
    strcpy(configPath, argv[1]);
  }
  else {
    cout << "One argument required: YAML config file." << endl;
    exit(1);
  }
}

void CommandLineParser::print() {
  cout << "Command Line Arguments: " << endl;
  cout << "configPath = " << configPath << endl;
  options.printOptions();
}
