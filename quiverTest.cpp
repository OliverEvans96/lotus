#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Quiver.h"
#include "TCanvas.h"

using namespace std;

int main()
{
	ifstream f1("step4.txt");
	ifstream f2("step4.txt");

	int N=0;
	string s;
	
	while(!f1.eof())
	{
		getline(f1,s);
		N++;
	}

	//cout << "N=" << N << endl;

	double coords[3];
	double velocities[3];

	double *x = new double[N];
	double *y = new double[N];
	double *z = new double[N];

	double *vx = new double[N];
	double *vy = new double[N];
	double *vz = new double[N];

	double *r = new double[N];
	double *vr = new double[N];

	int k = 0;
	
	TCanvas *c = new TCanvas();
	Quiver *q = new Quiver(10,0,150,10,25,50);
	q->Reset();

	//cout << "Read" << endl;
	while ((!f2.eof()) && (true))
	{
		f2 >> x[k] >> y[k] >> z[k] >> vx[k] >> vy[k] >> vz[k];
		f2.ignore(256,'\n');

		r[k] = sqrt(x[k]*x[k]+y[k]*y[k]);
		if(r[k] != 0)
			vr[k] = 2/r[k]*(x[k]*vx[k]+y[k]*vy[k]);
		else
			vr[k] = 0;

		q->Fill(r[k],z[k],vr[k],vz[k]);

		//cout << r[k] << " " << vr[k] << " " << z[k] << " " << vz[k] << endl;
		k++;
	}

	//cout << "Q" << endl;


	//cout << "Title" << endl;
	q->SetTitle("Plot!");
	//cout << "Draw" << endl;
	q->SetArrowParams(40,0.03,2,0.5);
	q->Draw(c);
	//cout << "Save" << endl;
	q->SaveAs("quiver.png");

	delete [] x,y,z,r,vx,vy,vz,vr;
	delete q,c;

	return 0;
}
