void horiz()
{
	cout << "Loaded!" << endl;
	return 0;
}

TGraph *horizontalHist(TH1D* hist)
{
	char *title, *xLabel, *yLabel;
	int n;
	double *x,*y;

	//Copy points
	TGraph* g =	new TGraph(hist);
	n=g->GetN();
	x=g->GetX();
	y=g->GetY();

	//Copy title and axis labels
	title=(char*) hist->GetTitle();
	xLabel=(char*) hist->GetXaxis()->GetTitle();
	yLabel=(char*) hist->GetYaxis()->GetTitle();

	//Flip x & y
	TGraph* g1 = new TGraph(n,y,x);

	//Copy titles
	g1->SetTitle(title);
	g1->GetXaxis()->SetTitle(yLabel);
	g1->GetYaxis()->SetTitle(xLabel);

	//Copy line attributes
	g1->SetLineStyle(g->GetLineStyle());
	g1->SetLineWidth(g->GetLineWidth());
	g1->SetLineColor(g->GetLineColor());

	//Copy fill attributes
	g1->SetFillStyle(g->GetFillStyle());
	g1->SetFillColor(g->GetFillColor());

	//Draw
	g1->Draw();

	//Delete
	//delete g;
	return g1;
}
