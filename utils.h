#ifndef DROPLET_ANALYSIS_UTILS_H
#define DROPLET_ANALYSIS_UTILS_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <sys/stat.h>
#include <algorithm>

//Split a string into a string vector of words
vector<string> strSplit(string str);

//Get index and coordinates from string array of words
void strToData(double *coords,double *velocities,double& dipole,string line);

//Choose the highest value from a vector
double max(vector<double> v);

//Find the mean of a double vector
double mean(vector<double> v);

//Square
double square(double x);

//Arctanh
double atanh(double x);

//Is a<b?
bool isLess(int a,int b);

//Is x in v?
bool isIn(int x, vector<int>v);

//Check whether a file exists
bool file_exists(const string& name);

//round up to nearest multiple of a number
int roundUp(int numToRound, int multiple);

//Find the minimum value of a vector or double pointer
double findMinimum(vector<double> v);

//Find the maximum value of a vector or double pointer
double findMaximum(vector<double> v);

#endif
