#ifndef PARSE_H
#define PARSE_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>

#include "Utils.h"

using namespace std;

////////// sortWaters

//Count the number of water atoms in the system
int countAtoms(ifstream &inFile);

//Count the number of timesteps
int countSteps(ifstream &inFile);

//Read data for one timestep from file
void readData(ifstream &inFile,string *header,int *atomNumbers,string *wholeLines,int numAtoms);

//Quick sort!
void quickSort(int n[],string s[],int first,int last);

//Write sorted data to new file
void writeData(ofstream &outFile,string *header,string *wholeLines,int numAtoms);

//Swap elements i and j of parallel arrays
void swap(int n[],string s[],int i,int j);

//Get the first word from a string as an integer
int getFirstInt(string str);

/////// substrateDensity

//Read data for one timestep
void readStep(ifstream &inFile,double *bounds,double *coords,int *elements,int numAtoms,int numSteps,int &timestep,string *arr);

//Calculate density(z) for one timestep
void calcDensity(int numAtoms,double *bounds,double *coords,int *elements,int nZ,double zWidth,double *zVals,double *density,double *center);

//Calculate the center of mass of the substrate for one timestep
void centerOfMass(int numAtoms,double *coords,int *elements,double *center);

//Write data for one timestep
void writeStep(ofstream &densFile,ofstream &centerFile,int nZ,double *density,double *center,int timestep);

//Split a string into a string array of words
string *strSplit(string str,string *arr);

//Get index and coordinates from string array of words
void strToData(double *&coords,int *&elements,string line,double bounds[6],string *arr);


//////// calculations_centered

//Read data for one timestep
void readStep(ifstream &inFile,int *OIndices,double *OCoords,double *prevOCoords,double *nextOCoords,int *HIndices,double *HCoords,double *nextHCoords,int numMolecules,int numSteps,int &timestep);

//Calculate velocities for one timestep
void calculateStepVelocities(double *nextOCoords,double *prevOCoords,double *velocities,double *speeds,double dt,int numMolecules);

//Calculate dipoles for one timestep
void calculateStepDipoles(double *OCoords,double *HCoords,double *dipoles,int numMolecules);

//Write data for one timestep
void writeStep(ofstream &outFile,int *OIndices,double *OCoords,double *velocities,double *speeds,double *dipoles,int numMolecules,int timestep);

//Get index and coordinates from string array of words
void strToData(int *&index,double *&coords,string line,double bounds[6]);

#endif
