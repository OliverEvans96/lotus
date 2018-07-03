/**
   @file Utils.h

   Generally useful functions and constants.
   These are not specific to one usage,
   and may be used in multiple files.
*/

#ifndef DROPLET_ANALYSIS_UTILS_H
#define DROPLET_ANALYSIS_UTILS_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cmath>
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
/// Multiply amu/AA^3 by this number to get g/cc.
const double NANO_DENS_TO_MACRO = 1.0/0.60221409;

/// Count the number of lines in a file.
int countLines(ifstream &inFile);

/// Elementwise vector addition.
vector<double> add(vector<double> u, vector<double> v);
/// Elementwise vector multiplication.
vector<double> mult(vector<double> u, vector<double> v);
/// Sum the elements in a vector.
double sum(vector<double> v);
/// Sample standard deviation of a vector. DOF = n-1.
double stddev(vector<double> v);
/// Maximum element of a vector.
double max(vector<double> v);
/// Mean of a vector.
double mean(vector<double> v);

/// Minimum element of an array of length @p n.
double min(double *v, int n);
/// Maximum element of an array of length @p n.
double max(double *v, int n);
/// Mean of an array of length @p n.
double mean(double *v, int n);
/// Sample standard deviation of an array of length @p n.
double stddev(double *v, int n);

/// Alias for @p x * @p x.
double square(double x);

/// Hyperbolic arctangent.
double atanh(double x);

/// Minimum of two numbers, @p a and @p b.
double min(double a, double b);
/// Maximum of two numbers, @p a and @p b.
double max(double a, double b);

/// Is @p a < @p b?
bool isLess(int a,int b);

/// Is @p x in @p v?
bool isIn(int x, vector<int>v);

/// Is @p substr a substring of @p str?
bool isIn(string str, string substr);

/// Check whether a file exists and is not a directory.
bool file_exists(const char* pathname);
/// Check whether a directory exists.
bool dir_exists(const char* pathname);

/// Round up to nearest multiple of a number
int roundUp(int numToRound, int multiple);

/** Given a TH1D and a two bin numbers,
    draw a line between the points and solve for
    the @f$x @f$ value where @f$y=y_c @f$
    (@p yc stands for ``y_cutoff``).
*/
double solveLinear(TH1D *hist,int bin1,int bin2,double yc);

/// Remove the trailing forward slash from @p path if it is there.
void stripTrailingSlash(char* strippedPath, const char* path);
/// Join two paths by a forward slash (using #stripTrailingSlash)
void joinPath(char* path, const char* prefix, const char* suffix);

#endif
