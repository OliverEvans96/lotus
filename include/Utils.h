/**
   @file Utils.h

   @brief Generally useful functions and constants.

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
/// @brief Multiply @f$\mathrm{amu}/\mathrm{\mathring{A}}^3@f$ by this number to get @f$\mathrm{g}/\mathrm{cm}^3 @f$.
const double NANO_DENS_TO_MACRO = 1.0/0.60221409;

/// @brief Count the number of lines in a file.
int countLines(ifstream &inFile);

/// @brief Elementwise vector addition.
vector<double> add(vector<double> u, vector<double> v);
/// @brief Elementwise vector multiplication.
vector<double> mult(vector<double> u, vector<double> v);
/// @brief Sum the elements in a vector.
double sum(vector<double> v);
/// @brief Sample standard deviation of a vector. DOF = n-1.
double stddev(vector<double> v);
/// @brief Maximum element of a vector.
double max(vector<double> v);
/// @brief Mean of a vector.
double mean(vector<double> v);

/// @brief Minimum element of an array of length @p n.
double min(double *v, int n);
/// @brief Maximum element of an array of length @p n.
double max(double *v, int n);
/// @brief Mean of an array of length @p n.
double mean(double *v, int n);
/// @brief Sample standard deviation of an array of length @p n.
double stddev(double *v, int n);

/// @brief Alias for @p x * @p x.
double square(double x);

/// @brief Hyperbolic arctangent.
double atanh(double x);

/// @brief Minimum of two numbers, @p a and @p b.
double min(double a, double b);
/// @brief Maximum of two numbers, @p a and @p b.
double max(double a, double b);

/// @brief Is @p a < @p b?
bool isLess(int a,int b);

/// @brief Is @p x in @p v?
bool isIn(int x, vector<int>v);

/// @brief Is @p substr a substring of @p str?
bool isIn(string str, string substr);

/// @brief Check whether a file exists and is not a directory.
bool file_exists(const char* pathname);
/// @brief Check whether a directory exists.
bool dir_exists(const char* pathname);

/// @brief Round up to nearest multiple of a number
int roundUp(int numToRound, int multiple);

/** @brief Given a TH1D and a two bin numbers,
    draw a line between the points and solve for
    the @f$x @f$ value where @f$y=y_c @f$
    (@p yc stands for ``y_cutoff``).
*/
double solveLinear(TH1D *hist,int bin1,int bin2,double yc);

/// @brief Remove the trailing forward slash from @p path if it is there.
void stripTrailingSlash(char* strippedPath, const char* path);
/// @brief Join two paths by a forward slash (using #stripTrailingSlash)
void joinPath(char* path, const char* prefix, const char* suffix);

#endif
