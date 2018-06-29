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
	~CircleFit();

  void setContext(SimData &simData, TGraph* _gCirclePoints);
  void setGraph(TGraph* _gCirclePoints);
  void updatePoints();
	double GetChi2s();
	void GuessFit();
	void Fit();
	double LinearResidual(const double *X);
	double SumOfSquares(const double *X);
	double Intersect(double c);
	double GetContactAngle();
	double LinearContactAngle();
	double GetHeight();
	void SetTangentLine(TLine *tangentLine);
	double getXCenter();
	double getYCenter();
	double getRadius();
	void Print();
  int GetNumPoints();

  void DeletePoints(vector<int> indices);
  double GetResidual(int i);
  void Refine();

	bool intersected; //Whether Intersect() has been called
  double gx0, gy0, gr; // Guessed points (from MLS)

private:
  SimData* simDataPtr;
  Options options;
  char fitOptions[16];

  TGraph *gCirclePoints;
	//Variables
  vector<double> x, y;
  int n; // number of points
	double x0,y0,r; //x-center,y-center, and radius
	double x1,y1; //Intersection of bulk and monolayer
	double cosTheta; //Cos(contact angle)
  double thetaDeg; // contact angle (degrees)
	double height; //Droplet height at x=0
	double m,b; //Parameters for tangent line
	double x2,y2,x3,y3; //Tangent line points

	//Fitting variables
	TMinuitMinimizer minimizer, linMin;
	double A,B,C,D,E;
	double sumsq;
	double width;
	double stepVal; //Step size value
	double step[3]; //Step size vector
	double init[3]; //Initial values
	double chi2s; //Error in circle fitting

	// Linear fit variables
	double cutoff;
	vector<double> xLinFit, yLinFit;
	double m_lin, b_lin;

	//Functions
	void findRadius();
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
  double solveLinear(int bin1, int bin2, double yc);
  void guessTanhFit();
  void initialGuess(double _ld=2.0, double _w=20.0, double _x0=50.0);
  bool isEmpty();
  void solve();
  double residual();
  bool good();

  double getBoundary();
  double getWidth();
  double getLiquidDensity();

  char rowOrCol[4];
  int rowColNum;

  SimData *simDataPtr;
  TF1 *fTanh;
  TH1D *hTanh;
};

#endif
