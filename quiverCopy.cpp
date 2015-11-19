// Vector field plot from two 2D histograms

#include <iostream>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>

#include "TStyle.h"
#include "TH2D.h"
#include "TH1I.h"
#include "TArrow.h"
#include "TCanvas.h"
#include "TText.h"

using namespace std;

double min(double x, double y)
{
	if (x<y)
		return x;
	else
		return y;
}

double minSet(double *x,int N)
{
	double min;

	//Don't want to return 0 unless all values are 0.
	for(int i=0;i<N;i++)
	{
		if( x[i] != 0)
		{
			min = x[i];
			break;
		}
	}


	for (int i=1;i<N;i++)
	{
		if ((x[i]<min) && (x[i] != 0))
		{
			min=x[i];
		}
	}
	
	return min;
}

double maxSet(double *x,int N)
{
	double max;

	// Don't want to return 0 unless all values are 0.
	for(int i=0;i<N;i++)
	{
		if(x[i] != 0)
		{
			max = x[i];
			break;
		}
	}

	for (int i=1;i<N;i++)
	{
		if (x[i]>max)
			max=x[i];
	}
	
	return max;
}

int main(int argc,char* argv[])
{

	TCanvas *c1 = new TCanvas("c1","c1",1000,800);

	const int N = 10;

	const int xlo = 1;
	const int xhi=4;
	const int nx = 3;

	const int ylo = 1;
	const int yhi=4;
	const int ny = 3;

	TH2D *h = new TH2D("h","h",nx,xlo,xhi,ny,ylo,yhi); // Dummy for plot
	TH2I *hCount = new TH2I("hCount","hCount",nx,xlo,xhi,ny,ylo,yhi); // Count entries
	TH2D *hvx = new TH2D("hvx","hvx",nx,xlo,xhi,ny,ylo,yhi); // v_x
	TH2D *hvy = new TH2D("hvy","hvy",nx,xlo,xhi,ny,ylo,yhi); // v_y

	double x[N] = {1,2,3,2,3,2,1,3,2,3};
	double y[N] = {3,2,1,3,2,3,2,3,1,3};
	double vx[N] = {1,-3,2,3,-2,3,-2,1,3,-2};
	double vy[N] = {1,-2,1,-3,2,-1,2,-3,2,1};

	// Tails
	double x1[N];
	double y1[N];

	// Velocities
	double vxk[N];
	double vyk[N];
	
	int nElements;
	double norm[N];
	double minWidth;
	double padding = 0.05; // % Padding on smallest box dimension
	double arrowLen;
	double factor;

	TArrow *a[nx*ny];
	TText *l[nx*ny];
	
	// Fill Histograms
	for (int i=0;i<N;i++)
	{
		hCount->Fill(x[i],y[i]);
		hvx->Fill(x[i],y[i],vx[i]);
		hvy->Fill(x[i],y[i],vy[i]);
	}

	// Calculate norms
	int k = 0;
	for (int i=0;i<nx;i++)
	{
		for (int j=0;j<ny;j++)
		{
			// Sums -> Averages (Divide by number of elements per bin)
			nElements = int(hCount->GetBinContent(i+1,j+1));

			if (nElements != 0)
			{
				vxk[k] = hvx->GetBinContent(i+1,j+1) / nElements;
				vyk[k] = hvy->GetBinContent(i+1,j+1) / nElements;
			}
			else
				vxk[k] = vyk[k] = 0;

			x1[k] = h->GetXaxis()->GetBinCenter(i+1);
			y1[k] = h->GetYaxis()->GetBinCenter(j+1);

			cout << "(" << x1[k] << "," << y1[k] << "): " << nElements << ", " << vxk[k] << " " << vyk[k] << endl;

			norm[k] = sqrt(vxk[k]*vxk[k]+vyk[k]*vyk[k]);

			k++;
		}
	}

	double minVal=minSet(norm,k);
	double maxVal=maxSet(norm,k);
	
	for (int i=0; i<k; i++)
	{
		cout << norm[i] << endl;
	}

	cout << "minMax " << minVal << " " << maxVal << endl;

	// Draw hist & color scale
	const int nLevels = 49;
	double levels[nLevels];
	int color;

	for (int i=1;i<nLevels;i++)
		levels[i] = minVal + (maxVal - minVal) / (nLevels - 1) * (i);
	levels[0]=0.01;
	cout << "sizeof(levels) = " << sizeof(levels) << endl;
	
	gStyle->SetOptStat(0);
	h->SetContour((sizeof(levels)/sizeof(double)), levels);
	h->SetBinContent(0,0,0);
	h->DrawCopy("col text"); // Draw axes
	h->GetZaxis()->SetRangeUser(minVal, maxVal);
	h->DrawCopy("z same"); // Draw color scale

	stringstream normStream;
	string normString;

	// Draw arrows
	k = 0;
	for (int i=0;i<nx;i++)
	{
		for (int j=0;j<ny;j++)
		{
			minWidth = min(h->GetXaxis()->GetBinWidth(i+1),h->GetYaxis()->GetBinWidth(j+1));
			arrowLen = (1.0-padding)*minWidth/2;
			factor = arrowLen/(norm[k]);

			if (norm[k] != 0)
			{
				a[k] = new TArrow(x1[k],y1[k],x1[k]+factor*vxk[k],y1[k]+factor*vyk[k],0.05,"|>");
				a[k]->SetAngle(40);
				a[k]->SetLineWidth(2);

				color = int((norm[k]-minVal)*nLevels/(maxVal-minVal));
				a[k]->SetFillColor(51+color);
				a[k]->SetLineColor(51+color);
				a[k]->Draw();

				//Convert norm to string for printing
				normStream << norm[k];
				normString = normStream.str();
				normStream.str("");

				l[k] = new TText(x1[k],y1[k],normString.data());
				l[k]->Draw();
			}

			k++;
		}

	}


	c1->SaveAs("quiver.png");

	return 0;
}

