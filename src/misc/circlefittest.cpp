#include <iostream>
#include <fstream>
#include "CircleFitClass.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TLegend.h"

using namespace std;

const int N_POINTS = 55;

int main()
{
	ifstream inFile("init_points_40A_half.txt");

	double x[N_POINTS];
	double y[N_POINTS];
	for(int i=0; i<N_POINTS;i++)
	{
		inFile >> x[i] >> y[i];
		inFile.ignore(256,'\n');
	}

	vector<double> x_vect(x,x+N_POINTS);
	vector<double> y_vect(y,y+N_POINTS);

	TCanvas *cv = new TCanvas();
	TH2D *h = new TH2D("h","h",200,0,120,200,25,60);
	cout << endl;
	CircleFit c1(h);
	c1.Define("name",x_vect,y_vect);
	cout << endl;
	TGraph *g = new TGraph(N_POINTS,x,y);
	g->SetMarkerStyle(20);

	g->Draw("ap");
	cout << "Draw" << endl;
	TEllipse *e1 = c1.Draw();
	e1->SetLineColor(kRed);

	g->Draw("same p");

	cv->SaveAs("CircleDataTest.png");
	cv->SaveAs("CircleDataTest.C");

	delete cv,g,e1;
}
