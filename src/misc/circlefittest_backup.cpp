#include <iostream>
#include <fstream>
#include "CircleFitClass.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TLegend.h"

using namespace std;

const int N_POINTS = 110;

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

	TCanvas *cv = new TCanvas();
	TH2D *h = new TH2D("h","h",200,0,120,200,25,60);
	cout << endl;
	CircleFit c1("CircleData_frame1.txt",h,"unif");
	cout << endl;
	CircleFit c2("CircleData_frame1.txt",h,"kasa");
	cout << endl;
	CircleFit c3("CircleData_frame1.txt",h,"mls");
	cout << endl;
	TGraph *g = new TGraph(N_POINTS,x,y);
	g->SetMarkerStyle(20);

	g->Draw("ap");
	TEllipse *e1 = c1.Draw();
	e1->SetLineColor(kRed);
	TEllipse *e2 = c2.Draw();
	e2->SetLineColor(kBlue);
	TEllipse *e3 = c3.Draw();
	e3->SetLineColor(kGreen);

	g->Draw("same p");

	TLegend *l = new TLegend(.65,.7,.85,.85);
	l->AddEntry(e1,"unif","l");
	l->AddEntry(e2,"kasa","l");
	l->AddEntry(e3,"mls","l");
	l->Draw();

	cv->SaveAs("CircleDataTest.png");
	cv->SaveAs("CircleDataTest.C");
	delete cv,g,e1,e2,e3;
}
