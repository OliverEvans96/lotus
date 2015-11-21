#ifndef QUIVER_H
#define QUIVER_H

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>

#include "TStyle.h"
#include "TH2D.h"
#include "TH1I.h"
#include "TArrow.h"
#include "TCanvas.h"
#include "TText.h"
#include "TMarker.h"

using namespace std;

class Quiver
{
public:
	Quiver(int nx, double xlo, double xhi, int ny, double ylo, double yhi);
	~Quiver();
	void Fill(double x, double y, double vx, double vy);
	void Calculate();
	void Draw(TCanvas *c=gPad->GetCanvas(),int pad=gPad->GetNumber());

	void Reset();
	void SetLevels(double min, double max);
	void SetDrawOptions(bool drawMarkers,bool drawText);
	void SetArrowParams(double arrowAngle, double arrowHead, int arrowWidth, double padding);
	void SetTitle(const char *title);
	void SaveAs(const char *filename);

private:
	// Functions
	double min(double x, double y);
	double minSet(vector<double> x, int N);
	double maxSet(vector<double> x, int N);

	// Variables

	// Basic - Data
	int N;
	int k;

	int nx;
	double xlo;
	double xhi;

	int ny;
	double  ylo;
	double yhi;

	// Canvas & Hists
	bool canvasCreated;
	TCanvas *cQuiver;
    TH2D *h;
	TH2I *hCount;
	TH2D *hvx;
	TH2D *hvy;

	// Calculations
	
	// Positions
	vector<double> x1;
	vector<double> y1;

	// Velocities
	vector<double> vxk;
	vector<double> vyk;

	// Others
	int nElements;
	vector<double> norm;
	double minWidth;
	double arrowLen;
	double factor;
	double minVal;
	double maxVal;

	// Draw  options
	bool drawMarkers;
	bool drawText;

	// Arrow Parameters
	double arrowAngle;
	double arrowHead;
	int arrowWidth;
	double padding;

	// Levels & color
	bool setLevels;
	int nLevels;
	double *levels;
	int color;

	// String
	stringstream normStream;
	string normString;

	// Arrows & Text
	TArrow **a;
	TText **l;
	TMarker **m;

	double xLen;
	double yLen;
};


#endif
