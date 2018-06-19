#ifndef DROPLET_ANALYSIS_UTILS_H
#define DROPLET_ANALYSIS_UTILS_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <map>
#include <sstream>
#include <sys/stat.h>

#include "TCanvas.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TText.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TError.h"

using namespace std;

const double PI = 3.141592653589793;
// Multiply amu/AA^3 by this number to get g/cc
const double NANO_DENS_TO_MACRO = 1.0/0.60221409;

//Choose the highest value from a vector
double max(vector<double> v);

//Find the mean of a double vector
double mean(vector<double> v);
double mean(double *v, int n);

//Square
double square(double x);

//Arctanh
double atanh(double x);

//Is a<b?
bool isLess(int a,int b);

//Is x in v?
bool isIn(int x, vector<int>v);

// Does the string contain a substring?
bool isIn(string str, string substr);

//Check whether a file exists
bool file_exists(const string& name);

//round up to nearest multiple of a number
int roundUp(int numToRound, int multiple);

//Find the minimum value of a vector or double pointer
double findMinimum(vector<double> v);

//Find the maximum value of a vector or double pointer
double findMaximum(vector<double> v);

//Given a TH1D and a two bin numbers, draw a line between the points and solve for where y=yc (y_cutoff)
double solveLinear(TH1D *hist,int bin1,int bin2,double yc);

#endif
