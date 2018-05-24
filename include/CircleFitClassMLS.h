#ifndef CIRCLEFITCLASSMLS_H
#define CIRCLEFITCLASSMLS_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include "TCanvas.h"
#include "TEllipse.h"
#include "TGraph.h"
#include "TEllipse.h"
#include "TLine.h"

using namespace std;

class CircleFit
{
public:
	CircleFit();
	CircleFit(char* filename);
	CircleFit(char* name,vector<double> xCoords,vector<double> yCoords);
	~CircleFit();
	
	void Define(char* name,vector<double> xCoords,vector<double> yCoords);
	void GetSize();
	void Fill(char* filename);
	void Fit();
	double Intersect(double c);
	double ContactAngle();
	double Height();
	TEllipse* Draw(bool drawPoints=false);
	TLine* DrawTangentLine();
	double GetXCenter();
	double GetYCenter();
	double GetRadius();
	void Print();
		
private:
	//Variables
	double PI; //Pi!
	char* name; //Circle name
	vector<double> x,y; //List of points
	int n; //Number of points
	double x0,y0,r; //x-center,y-center, and radius
	double x1,y1; //Intersection of bulk and monolayer
	double theta; //Contact angle
	double cosTheta; //Cos(contact angle)
	double height; //Droplet height at x=0
	bool intersected; //Whether Intersect() has been called
	double m,b; //Parameters for tangent line
	double margins[4]; //Margins of canvas
	double xlo,ylo,xhi,yhi; //Edges of histogram
	double x2,y2,x3,y3; //Tangent line points
	double minAngle,maxAngle; //Range over which to draw circle

	//Fitting variables
	double A,B,C,D,E;
	double tmpx,tmpy,tmpr;
	
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
