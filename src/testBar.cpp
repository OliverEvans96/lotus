#include <iostream>
#include <cmath>
#include <iostream>
#include "TStyle.h"
#include "TH2D.h"
#include "TArrow.h"
#include "TCanvas.h"


int main()
{
	const Double_t min = 0.9;
	const Double_t max = 1.1;

	const Int_t nLevels = 10;
	Double_t levels[nLevels];

	for(int i = 1; i < nLevels; i++) {
	levels[i] = min + (max - min) / (nLevels - 1) * (i);
	}
	levels[0] = 0.01;

	TCanvas *c1 = new TCanvas();
	TH2D *h  = new TH2D("h", "h", 10, 0, 10, 10, 0, 10);
	gStyle->SetOptStat(0);
	h->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
	h->SetBinContent(0,0,0);
	h->DrawClone("col text"); // Draw axes
	h->GetZaxis()->SetRangeUser(min, max);
	h->Draw("z same"); // Draw color scale

	c1->SaveAs("test.png");
}
