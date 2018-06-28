#ifndef FITTING_H
#define FITTING_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include "TCanvas.h"
#include "TEllipse.h"
#include "TFile.h"
#include "TGraph.h"
#include "TEllipse.h"
#include "TH2D.h"
#include "TLine.h"
#include "TMinuitMinimizer.h"
#include "Math/Functor.h"

#include "Utils.h"
#include "MDBase.h"
#include "Parameters.h"

using namespace std;

class CircleFit {
public:
  CircleFit();
	CircleFit(char* filename,TH2D *givenHist);
	CircleFit(char* name,vector<double> xCoords,vector<double> yCoords);
	~CircleFit();

  void setHist(TH2D *givenHist);
	void Define(char* name,vector<double> xCoords,vector<double> yCoords);
	double GetChi2s();
	void GetSize();
	void Fill(char* filename);
	void GuessFit();
	void Fit();
	double LinearResidual(const double *X);
	double SumOfSquares(const double *X);
	double Intersect(double c);
	double ContactAngle();
	double LinearContactAngle();
	double Height();
	void DrawBadPoints(TGraph *givenBadPointsGraph);
	void AddBadPoint(double xp,double yp);
	TEllipse* Draw(bool drawPoints=false);
	TLine* DrawTangentLine();
	double getXCenter();
	double getYCenter();
	double getRadius();
	void Print();
  int GetNumPoints();

	bool intersected; //Whether Intersect() has been called
  double gx0, gy0, gr; // Guessed points (from MLS)

private:
	//Variables
	char* name; //Circle name
	vector<double> x,y; //List of points
	int n; //Number of points
	double x0,y0,r; //x-center,y-center, and radius
	double x1,y1; //Intersection of bulk and monolayer
	double theta; //Contact angle
	double cosTheta; //Cos(contact angle)
  double thetaDeg; // contact angle (degrees)
	double height; //Droplet height at x=0
	double m,b; //Parameters for tangent line
	double margins[4]; //Margins of canvas
	double xlo,ylo,xhi,yhi; //Edges of histogram
	double x2,y2,x3,y3; //Tangent line points
	double minAngle,maxAngle; //Range over which to draw circle

	//Histogram variables
	TH2D *hist; //Associated histogram
	int xbin,ybin; //Bin containing point

	//Fitting variables
	TMinuitMinimizer minimizer, linMin;
	double A,B,C,D,E;
	double sumsq;
	double width;
	double stepVal; //Step size value
	double step[3]; //Step size vector
	double init[3]; //Initial values
	double tmpx,tmpy,tmpr;
	double chi2s; //Error in circle fitting

	// Linear fit variables
	double cutoff;
	vector<double> xLinFit, yLinFit;
	double m_lin, b_lin;


	//Functions
	double sum(vector<double> v);
	double square(double x);
	vector<double> add(vector<double> u,vector<double> v);
	vector<double> mult(vector<double> u,vector<double> v);
	void findRadius();
	void countLines(ifstream &inFile);
	double mean(vector<double> x);
	double stddev(vector<double> x);
	double atanh(double x);
	bool inGraph(TGraph *g,double xCheck,double yCheck);
};

class TanhFit {
  Options options;
  char fitOptions[16];

  //Set Bounds on parameters
  // TODO: Set bounds from options
  double fitBounds[6];

  double ld;
  double w;
  double x0;

  int err;
  bool empty;

public:
  TanhFit();
  ~TanhFit();

  void setContext(SimData &simData);
  void createFunction();
  void setHist(TH1D* _hTanh);
  void setFitBounds();
  void setFitType(const char* _rowOrCol);
  void setFitNum(int num);
  void initialGuess(double _ld=2.0, double _w=20.0, double _x0=50.0);
  bool isEmpty();
  void solve();
  double residual();
  bool good();

  double getBoundary();
  double getWidth();
  double getLiquidDensity();
  SimData *simDataPtr;

  char rowOrCol[4];
  int rowColNum;

  TF1 *fTanh;
  TH1D *hTanh;
};

#endif
