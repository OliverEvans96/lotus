#ifndef CIRCLEFITCLASS_H
#define CIRCLEFITCLASS_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include "TCanvas.h"
#include "TEllipse.h"
#include "TGraph.h"

using namespace std;

class CircleFit
{
public:
	CircleFit();
	CircleFit(char* filename);
	CircleFit(char* name,vector<double> xCoords,vector<double> yCoords);
	~CircleFit();
	
	void GetSize();
	void Fill(char* filename);
	void Fit();
	double Intersect(double c);
	double ContactAngle();
	double Height();
	void Draw(bool drawPoints=false);
	double GetXCenter();
	double GetYCenter();
	double GetRadius();
	void Print();
		
private:
	//Variables
	char* name;
	vector<double> x,y; //List of points
	int n; //Number of points
	double x0,y0,r; //x-center,y-center, and radius
	double x1,y1; //Intersection of bulk and monolayer
	double theta; //Contact angle
	double cosTheta; //Cos(contact angle)
	double height; //Droplet height at x=0
	bool intersected; //Whether Intersect() has been called
	
	//Functions
	double sum(vector<double> v);
	double square(double x);
	vector<double> add(vector<double> u,vector<double> v);
	vector<double> mult(vector<double> u,vector<double> v);
	void findRadius();
	void countLines(ifstream &inFile);
	double atanh(double x);
	

};

#endif
