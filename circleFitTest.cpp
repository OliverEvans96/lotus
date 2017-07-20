#include <iostream>
#include "CircleFitClass.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"

using namespace std;

int main()
{
	CircleFit c("circlepoints.txt");
	TGraph *g = new TGraph("circlepoints.txt");

	TCanvas *can = new TCanvas();
	g->GetXaxis()->SetRangeUser(-150,150);
	g->GetYaxis()->SetRangeUser(-150,150);
	g->SetMarkerStyle(20);
	g->SetMarkerSize(0.5);
	g->Draw("AP");
	c.Draw();
	c.Print();
	can->SaveAs("circletest.png");
	can->SaveAs("circletest.C");

	return 0;
}

