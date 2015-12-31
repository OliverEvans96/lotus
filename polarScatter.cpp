//Read each timestep and plot z(r)

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <sys/stat.h>
#include <algorithm>
#include "TCanvas.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include "TFile.h"

#include "CircleFitClass.h"
#include "Quiver/Quiver.h"

using namespace std;

const double PI = 3.141592653589793;

//Count the number of water atoms in the first timestep
int countAtoms(ifstream &inFile);

//Split a string into a string array of words
vector<string> strSplit(string str);

//Get index and coordinates from string array of words
void strToData(double *coords,double *velocities,double& dipole,string line);

//Choose the highest value from a vector
double max(vector<double> v);

//Find the mean of a double vector
double mean(vector<double> v);

//Square
double square(double x);

//Arctanh
double atanh(double x);

//Guess boundary of water molecule by counting for a single row
double guessRowBoundary(TH2D* hist,int j);

//Guess boundary of water molecule by counting for a single column
double guessColBoundary(TH2D* hist,int i);

//Draw a TH1D horizontally, returning a TGraph which should probably be deleted
TGraph *horizontalHist(TH1D* hist);

//Find boundary points by tanh fitting for each row
TGraph* findBoundaryPoints(TH2D* hist,char* aOrR,double& xBulkMax,double &xMonoEdge,int stepNum,int timestep);

//Given a graph of points, a maximum x value, and the y coordinate of the interface with the substrate, fit a circle to the points
CircleFit fitCircle(TGraph* g,double xMax,int stepNum,int timestep);

//Calculate MSD - mean squared displacement
vector<double> calculateMSD(vector<double> xi,vector<double> yi,vector<double> zi,vector<double> x, vector<double> y,vector<double> z);

//Count the number of atoms which have migrated from the bulk to the monolayer this timestep
double bulkMonoExchange(vector<double>z,double zInterface,int stepsPerFrame);

//Keep track of which atoms join the monolayer, and save their radius scaled by the base radius
void joinMonoPosition(vector<double> r,vector<double>z,double zInterface,double baseRadius,TH1D* rScaledJoin);

//Keep track of which atoms leave the monolayer, and save their radius scaled by the base radius to the vector rScaledLeave
void leaveMonoPosition(vector<double> r,vector<double>z,double zInterface,double baseRadius,TH1D* rScaledLeave);

//Is a<b?
bool isLess(int a,int b);

//Is x in v?
bool isIn(int x, vector<int>v);

//Plot bulk and mono radii over time
void plotRadii(int numSteps,double* bulkEdge,double* monoEdge,char* name);

//Find the minimum value of a vector or double pointer
double findMinimum(vector<double> v);

//Find the maximum value of a vector or double pointer
double findMaximum(vector<double> v);

//Count the number of timesteps
int countSteps(ifstream &inFile);

//Find the lowest atom ID for a water molecule in this simulation
int lowestID(ifstream &inFile);

//Check whether a file exists
bool file_exists(const string& name);

//Read each timestep and plot 3d scatter plot
int main(int argc,char* argv[])
{	
	//Command line arguments: inLoc,outLoc
	cout << "Received " << argc+1 << " command line arguments." << endl;

	char* inLoc=argv[1];
	char* outLoc=argv[2];
	
	//File Number
	string fileName(inLoc);
	string halfName=fileName.substr(fileName.find("atom")).substr(4);
	string fileNumberStr=halfName.substr(0,halfName.find("/"));
	stringstream fileNumberSS;
	fileNumberSS << fileNumberStr;
	int fileNumber;
	fileNumberSS >> fileNumber;
	cout << "Reading file number " << fileNumber << endl;

	cout << "Open Stream" << endl;
	//Input Stream
	ifstream inFile(inLoc);

	if (inFile.good())
		cout << "Successfully opened " << inLoc << endl;
	else
		cout << "Failed to open " << inLoc << endl;
	
	//Output file streams
	FILE* avgStepData = fopen("avgStepData.txt","w");
 	fprintf(avgStepData,"%7.7s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s\n","#timestep","bulkEdge","monoEdge","contactAngle","dropletHeight","avgRadius","avgDipole","monoExchange","MSDx","MSDy","MSDz","MSDavg");
	fflush(avgStepData);
	
	FILE* instStepData = fopen("instStepData.txt","w");
 	fprintf(instStepData,"%7.7s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s \n","#timestep","avgRadius","avgDipole","monoExchange","MSDx","MSDy","MSDz","MSDavg");
	fflush(instStepData);
	
	//Variables
	string line;
	string input;
	int atomNum=0;
	int lineNum=1;
	double coords[3];
	double velocities[3];
	double dipole;
	int timestep; //Real step number for simulation
	int stepNum=0; //Step number relative to this file starting with 0
	int frameNum=0;
	bool loopFlag1=true;
	bool loopFlag2=true;
	
	//Find out some things about the file
	cout << "Analyzing input file" << endl;
	ifstream testFile(inLoc);
	int numSteps=countSteps(testFile);
	int numAtoms=countAtoms(testFile);
	int firstID=lowestID(testFile);
	testFile.close();
	
	//Number of timesteps to average together to make a frame
	int stepsPerFrame=5;
	
	//Number of frames (collections of timesteps)
	//If numSteps%stepsPerFrame!=0, then the last frame will have fewer than the rest
	int numFrames=(int) ceil(numSteps/stepsPerFrame);

	//Coordinates
	vector<double> x(numAtoms);
	vector<double> y(numAtoms);
	vector<double> z(numAtoms);
	vector<double> r(numAtoms);
	vector<double> p(numAtoms);

	vector<double> vx(numAtoms);
	vector<double> vy(numAtoms);
	vector<double> vz(numAtoms);
	vector<double> vr(numAtoms);

	vector<double> cosT(numAtoms);
	
	//Initial positions
	vector<double> xi(numAtoms);
	vector<double> yi(numAtoms);
	vector<double> zi(numAtoms);
	
	//Read from file
	bool initWritten=file_exists("../init.txt");
	if (initWritten)
	{
		ifstream initIn("../init.txt");
		for(int i=0;i<numAtoms;i++)
		{
			initIn >> xi[i] >> yi[i] >> zi[i];
			initIn.ignore(256,'\n');
		}
		initIn.close();
	}
	
	//Center of mass
	double x0,y0;
	bool centroidKnown;
	
	//Step variables
	double xBulkMax;
	double xBulkEdge=10; //Start with a small value for joinMonoPosition()
	double xMonoEdge;
	double dropletHeight,contactAngle;
	double avgRadius,avgDipole;
	vector<double> msd;
	double exchange;
	double min;
	
	//Hist of where atoms join the monolayer
	TH1D* rScaledJoin = new TH1D("rScaledJoin","rScaledJoin",30,0,1.5);
	
	//Hist of where atoms leave the monolayer
	TH1D* rScaledLeave = new TH1D("rScaledLeave","rScaledLeave",30,0,1.5);
	
	double dz=1.5;
	double zlo=28.5;
	double zhi=160.5;
	int nz=(int) round((zhi-zlo)/dz);
	
	//Check that dz, zlo, and zhi are compatible
	if( !(fmod((zhi-zlo),dz)==0) )
	{
		cout << "(zhi-zlo) is not divisible by dz" << endl;
		return 1;
	}
	
	cout << "nz=" << nz << endl;
	
	double dA=500;
	double alo=0;
	double ahi=125500;
	int nA=(int) round((ahi-alo)/dA);
	
	double dr=dz;
	double rlo=sqrt(alo/PI);
	double rhi=sqrt(ahi/PI);
	int nr=(int) round((rhi-rlo)/dr);

	double vlo=0;
	double vhi=1e7;
	int nV=(int) round((vhi-vlo)/dA);

	double dp=dr;
	double plo=pow((3*vlo/(4*PI)),(1.0/3));
	double phi=pow((3*vhi/(4*PI)),(1.0/3));
	int  np=(int) round((phi-plo)/dp);
	
	//Bin volume
	double dV=dA*dz;
	
	//Check that dA, alo, and ahi are compatible
	if( !(fmod((ahi-alo),dA)==0) )
	{
		cout << "(Ahi-Alo) is not divisible by dA" << endl;
		return 1;
	}
	
	cout << "nA=" << nA << endl;
	
	//Conversion factors
	double gramsPerWaterMolecule=2.9914e-23;
	double A3perCC=1.0e24;
	double convFact=gramsPerWaterMolecule*A3perCC;

	//Skip first 2 lines
	for(int i=0;i<2;i++) inFile.ignore(256,'\n');
	lineNum+=2;
	
	//Fix canvas height & scale width to preserve aspect ratio
	int cH=1200;
	int cW=(int) round(cH*(rhi-rlo)/(zhi-zlo));
	TCanvas *cD = new TCanvas("cD","cD",cW,cH); //Dipole
	TCanvas *cA = new TCanvas("cA","cA",cW,cH); //2D density plot
	TCanvas *cQ = new TCanvas("cQ","cQ",cW,cH); //Quiver
	TCanvas *cVr = new TCanvas("cVr","cVr",cW,cH); //Plot v_r(z)
	TCanvas *cVp = new TCanvas("cVp","cVp",cW,cH); //Plot v_p(z)

	TCanvas *cAll = new TCanvas("cAll","cAll",cW,cH); //All 4 plots together
	cAll->Divide(2,2,0,0);

	//Values to use for equal-area bins
	double aVals[nA];
	double vVals[nA];
	double rVals[nA];
	double pVals[nA];
	for(int i=0;i<=nA;i++)
	{
		aVals[i]=i*dA;
		vVals[i]=i*dV;
		rVals[i]=sqrt(aVals[i]/PI);
		pVals[i]=pow(((3*vVals[i])/(4*PI)),1.0/3);
	}
		
	//Create Histograms
	TH2D *hA = new TH2D("hA","hA",nA,rVals,nz,zlo,zhi);
	TH1D *hD = new TH1D("hD","hD",20,0,1);
	TH1D *hVr = new TH1D("hVr","hVr",nz,zlo,zhi); //For tracking v_r(z)
	TH1D *hVp = new TH1D("hVp","hVp",nz,0,100); //For tracking v_rho(z) (rho measured from (0,0)
	TH1D *hZ = new TH1D("hZ","hZ",nz,zlo,zhi); //For counting # of atoms in each z bin
	TH1D *hP = new TH1D("hP","hP",np,zlo,zhi); //For counting # of atoms in each p bin

	//misc plot settings
	hA->SetStats(0);

	hVr->SetStats(0);
	hVr->SetLineWidth(2);
	hVr->GetXaxis()->SetTitle("z (#AA)");
	hVr->GetXaxis()->CenterTitle();
	hVr->GetYaxis()->SetTitle("v_{r} (#AA/fs)");
	hVr->GetYaxis()->CenterTitle();

	//Use OpenGL for antialiasing
	gStyle->SetCanvasPreferGL(true);
	
	//Quiver
	Quiver *q = new Quiver(nr/2,rlo,rhi,nz/2,zlo,zhi);
	q->SetArrowParams(30.0,0.01,2,0.05);
	q->SetLevels(0,2e-3);
	
	//Axis labels
	hA->GetXaxis()->SetTitle("r (#AA)");
	hA->GetXaxis()->CenterTitle();
	hA->GetYaxis()->SetTitle("z (#AA)");
	hA->GetYaxis()->CenterTitle();
	
	
	hD->GetXaxis()->SetTitle("cos#theta (rad)");
	hD->GetXaxis()->CenterTitle();
	hD->GetYaxis()->SetTitle("Occurrence");
	hD->GetYaxis()->CenterTitle();

	//Z boundary between monolayer and bulk (boundary between first and second bins
	double zInterface = zlo+dz;
	
	
	double colzMin=0.0;
	double colzMax=1.5;

	hA->SetMinimum(colzMin);
	hA->SetMaximum(colzMax);

// 	Edge r coordinates
	double bulkEdge[numFrames];
	double monoEdge[numFrames];
	
	//Set all values to zero
	for(int i=0;i<numSteps;i++)
	{
		bulkEdge[i]=0;
		monoEdge[i]=0;
	}
	
	//Determine whether to track mono atoms
	bool trackMonoAtoms=file_exists("../monoIDs.txt");
	int numMonoIDs=0;
	int id;
	vector<int> monoIDs;
	double monoR,monoZ;
	if(trackMonoAtoms)
	{
		cout << "Tracking allocating!" << endl;
		
		//Read IDs of atoms to track
		ifstream monoIDsIn("../monoIDs.txt");
		while(monoIDsIn.good())
		{
			monoIDsIn >> id;
			monoIDs.push_back((id-firstID)/3);
			numMonoIDs++;
			cout << "ID: " << id << " -> " << (id-firstID)/3 << endl;
		}
		monoIDsIn.close();
	}

	//Declare variables for velocity plot
	double sum,num;
	
	//Declaring graphs in the if block seems to produce errors
	TGraph* monoGraph1 = new TGraph(numMonoIDs); //Outline
	TGraph* monoGraph2 = new TGraph(numMonoIDs); //Fill
	
	if(trackMonoAtoms)
	{
		//Create graphs to draw to
		monoGraph1->SetMarkerStyle(20);
		monoGraph1->SetMarkerColor(kBlack);
		monoGraph1->SetMarkerSize(1.6);
		
		monoGraph2->SetMarkerStyle(20);
		monoGraph2->SetMarkerColor(kWhite);
		monoGraph2->SetMarkerSize(1.3);
	}
	
	//Create directories if they don't exist
	system("mkdir -p img/hist"); 
	system("mkdir -p img/mono");
	system("mkdir -p img/quiver");
	system("mkdir -p img/dipole");
	system("mkdir -p img/vr");
	system("mkdir -p img/vp");
	system("mkdir -p img/all");


	//Read data
	//For each timestep
	while(!inFile.eof()&&loopFlag1)
	{
		//Read the header
		inFile >> line >> timestep;
		inFile.ignore();
		lineNum++;
		
		cout << "Timestep " << timestep << " @ line " << lineNum-1 << endl;
		
		loopFlag2=true;
		
		atomNum=0;

		//Read atom data
		while(loopFlag2)
		{
			getline(inFile,line);
			lineNum++;
			if(line=="") loopFlag2=false;
			else
			{
				//Save data
				strToData(coords,velocities,dipole,line);
				x[atomNum]=coords[0];
				y[atomNum]=coords[1];
				z[atomNum]=coords[2];

				vx[atomNum]=velocities[0];
				vy[atomNum]=velocities[1];
				vz[atomNum]=velocities[2];

				cosT[atomNum]=dipole;
				atomNum++;
			}
		}
		
		//Save first timestep
		if( (!initWritten) && (stepNum==0) )
		{
			ofstream initOut("../init.txt");
			for(int i=0;i<numAtoms;i++)
			{
				xi[i]=x[i];
				yi[i]=y[i];
				zi[i]=z[i];
				
				//Save to file
				initOut << x[i] << " " << y[i] << " " << z[i] << endl;
			}
			initOut.close();
		}
		
		//Center
		//If the center has already been determined for this simulation, use that.
		//Otherwise, write it.
		if(!centroidKnown)
		{
			x0=mean(x);
			y0=mean(y);
			centroidKnown=true;
			
			if(file_exists("../centroid.txt"))
			{
				ifstream centroidIn("../centroid.txt");
				centroidIn >> x0 >> y0;
				centroidIn.close();
			}
			
			else
			{
				ofstream centroidOut("../centroid.txt");
				centroidOut << x0 << " " << y0;
				centroidOut.close();
			}
		
		}
		
		//Minimum
		min=findMinimum(z);

		//Track mono atoms
		if(trackMonoAtoms)
		{
			for(int i=0;i<numMonoIDs;i++)
			{
				monoR=sqrt(square(x[monoIDs[i]]-x0)+square(y[monoIDs[i]]-y0));
				monoZ=z[monoIDs[i]];
				//Fill graphs with radial positions of mono atoms
				monoGraph1->SetPoint(i,monoR,monoZ);
				monoGraph2->SetPoint(i,monoR,monoZ);
			}
		}
			
		//Fill Histograms
		for(int i=0;i<numAtoms;i++)
		{
			r[i]=sqrt(square(x[i]-x0)+square(y[i]-y0));
			p[i]=sqrt(square(r[i])+square(z[i]-zlo));

			vr[i]=(x[i]*vx[i]+y[i]*vy[i])/r[i]; // Chain rule

			hA->Fill(r[i],z[i],convFact/(dV*stepsPerFrame));
			hD->Fill(dipole,1/numAtoms);
			hVr->Fill(z[i],vr[i]);
			hVp->Fill(p[i],vr[i]);
			hZ->Fill(z[i]);
			q->Fill(r[i],z[i],vr[i],vz[i]);
			//cout << "Filling q: (" << vr[i] << "," << vz[i] << ") @ (" << r[i] << "," << z[i] << ")" << endl;

		}
//		cout << "min p = " << findMinimum(p) << endl;
//		cout << "max p = " << findMaximum(p) << endl;

		
// 		//Polar Scatter Plot
// 		ofstream psOut("pS.txt");
// 		for(int i=0;i<numAtoms;i++)
// 		{
// 			psOut << r[i] << " " << z[i] << endl;
// 			psOut.flush();
// 		}
// 		psOut.close();
		
		
		//Calculate average radius
		avgRadius=mean(r);
		
		//Calculate average dipole moment angle
		avgDipole=mean(cosT);
		
		//Calculate MSD
		msd=calculateMSD(xi,yi,zi,x,y,z);
		
		fprintf(instStepData,"%08d %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f  %15.6f\n",timestep,avgRadius,avgDipole,exchange,msd[0],msd[1],msd[2],msd[3]);
		
		//Update list of where atoms joined and left the monolayer
		joinMonoPosition(r,z,zInterface,xBulkEdge,rScaledJoin);
		leaveMonoPosition(r,z,zInterface,xBulkEdge,rScaledLeave);
		
		//Quiver
		stringstream title;
		title << "Radial Density (g/cc):  " << timestep;
		hA->SetTitle(title.str().data());
		q->SetTitle(title.str().data());

		
		//Calculations - Every frame (5 timesteps)
		if( (stepNum!=0) && ((stepNum+1)%stepsPerFrame==0) )
		{
			//Title
			stringstream title;
			title << "Radial Density (g/cc):  " << timestep;
			hA->SetTitle(title.str().data());

			title.str("");
			title << "Velocity Field: " << timestep;
			q->SetTitle(title.str().data());

			title.str("");
			title << "Dipole Moment Distribution: " << timestep;
			hD->SetTitle(title.str().data());

			title.str("");
			title << "2D Radial Velocity: " << timestep;
			hVr->SetTitle(title.str().data());

			title.str("");
			title << "3D Radial Velocity: " << timestep;
			hVp->SetTitle(title.str().data());

			//Draw dipole Histogram
			cD->cd();
			hD->Draw();
			title.str("");
			title << "img/dipole/step" << setw(8) << setfill('0') << timestep << ".png";
			cD->SaveAs(title.str().data());
			
			//Quiver
			title.str("");
			title << "img/quiver/step" << setw(8) << setfill('0') << timestep << ".png";
			q->Draw(cQ);
			cout << "Saving quiver" << endl;
			q->SaveAs(title.str().data());

			//2D Velocity Plot
			//Scale values to create average
			for(int i=0;i<nz;i++)
			{
				sum=hVr->GetBinContent(i+1);
				num=hZ->GetBinContent(i+1);
				if(num!=0)
					hVr->SetBinContent(i+1,abs(sum/num));
			//	cout << "hVr[" << i+1 << "] = " << hVr->GetBinContent(i+1) << endl;
			}
			//Plot
			cVr->cd();
			TGraph *g = horizontalHist(hVr);
			g->GetXaxis()->SetLimits(1e-6,1e-3);
			cVr->SetLogy();
			g->Draw("AL");
			cVr->Update();
			cVp->cd();
			//g->Draw("AL");
			q->Draw(cVr);
			title.str("");
			title << "img/vr/step" << setw(8) << setfill('0') << timestep << ".png";
			cVr->SaveAs(title.str().data());

			//3D Velocity Plot
			//Scale values to create average
			for(int i=0;i<nz;i++)
			{
				sum=hVp->GetBinContent(i+1);
				num=hP->GetBinContent(i+1); //NOT SURE HOW TO HANDLE THIS. hVp is kind of meaningless.
				if(num!=0)
					hVp->SetBinContent(i+1,abs(sum/num));
			}
			//Plot
			cVp->cd();
			TGraph *g1 = new TGraph(hVp);
			//hVr->Draw("AL");
			//g1->Draw("AL");
			g1->SetLineWidth(2);
			g1->GetXaxis()->SetTitle("#rho (#AA)");
			//g1->GetXaxis()->CenterTitle();
			g1->GetYaxis()->SetTitle("v_{r} (#AA/fs)");
			//g1->GetYaxis()->CenterTitle();

			title.str("");
			title << "img/vp/step" << setw(8) << setfill('0') << timestep << ".png";
			cVp->SaveAs(title.str().data());

			//Draw density histogram
			cA->cd();
			hA->Draw("colz");
			//Find boundary points
			TGraph* b2 = findBoundaryPoints(hA,"a",xBulkMax,xMonoEdge,stepNum,timestep);
			b2->SetMarkerStyle(20);
			b2->Draw("same L");
			stringstream aName;
			aName << "img/hist/abin" << setw(8) << setfill('0') << timestep << ".png";

			//Find interface & edge
			cout << endl;
			CircleFit Circle = fitCircle(b2,xBulkMax,stepNum,timestep);
			xBulkEdge=Circle.Intersect(zInterface);
			cout << "Bulk edge: " << xBulkEdge << endl;
			bulkEdge[frameNum]=xBulkEdge;
			cout << "Mono edge: " << xMonoEdge << endl;
			monoEdge[frameNum]=xMonoEdge;
			cA->SaveAs(aName.str().data());
			 	
			//Draw mono atoms
			if(trackMonoAtoms)
			{
				monoGraph1->Draw("P");
				monoGraph2->Draw("P");
				
				stringstream mName;
				mName << "img/mono/monoHist" << setw(8) << setfill('0') << timestep << ".png";
				cA->SaveAs(mName.str().data());
			}
		
			//Draw all 4 together
			cAll->cd(1);
			gPad->Clear();
			cA->DrawClonePad();
			cAll->cd(2);
			gPad->Clear();
			cVr->DrawClonePad();
			cAll->cd(3);
			gPad->Clear();
			cVp->DrawClonePad();
			cAll->cd(4);
			gPad->Clear();
			cQ->DrawClonePad();

			//Delete vr graph
			delete g;
			delete g1;

			title.str("");
			title << "img/all/step" << setw(8) << setfill('0') << timestep << ".png";
			cAll->SaveAs(title.str().data());
			
			//Find contact angle
			cout << "Intersect with " << zInterface << endl;
			contactAngle=Circle.ContactAngle();
			
			//Find droplet height by intersecting circle with y-axis
			dropletHeight=Circle.Height();
			
			//Calculate exchange
			exchange=bulkMonoExchange(z,zInterface,stepsPerFrame);
			
			//Find boundary
	// 		int numBoundaryPoints=b1->GetN();
	// 		double bPx[numBoundaryPoints],bPy[numBoundaryPoints];
	// 		
	// 		Check points:
	// 		cout << endl << "Points from main() b1" << endl;
	// 		for(int i=0;i<numBoundaryPoints;i++)
	// 		{
	// 			b1->GetPoint(i,bPx[i],bPy[i]);
	// 			cout << "("<<bPx[i]<<","<<bPy[i]<<")"<<endl;
	// 		}
	// 	

			hA->Reset();
			hD->Reset();
			q->Reset();
			
			//Save radii to file
			fprintf(avgStepData,"%08d %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n",timestep,xBulkEdge,xMonoEdge,contactAngle,dropletHeight,avgRadius,avgDipole,exchange,msd[0],msd[1],msd[2],msd[3]);
			fflush(avgStepData);
			
			frameNum++;
			
		}
		
		//Otherwise, fill exchange without expecting a result
		else
		{
			bulkMonoExchange(z,zInterface,stepsPerFrame);
		}
		
		//Only first timestep
		//loopFlag1=false;
		
		//Increment relative timestep counter
		stepNum++;
		
	}
	
	inFile.close();
	fclose(avgStepData);
	fclose(instStepData);
	
	
	//Join & Leave
	double I;
	TH1D* rExchange = rScaledJoin->Clone();
	//Normalize histograms
	rExchange->Add(rScaledLeave,-1);
	I=rScaledJoin->Integral();
	*rScaledJoin->Scale(1/I);
	I=rScaledLeave->Integral();
	*rScaledLeave->Scale(1/I);
	I=rExchange->Integral();
	*rExchange->Scale(1/I);
	TCanvas* cE = new TCanvas();
	cE->cd();
	rScaledJoin->Draw("");
	rScaledLeave->Draw("same");
	rExchange->Draw("same");
	cE->SaveAs("leaveJoinHists.png");
	
	ofstream leaveJoinFile("leavejoin.txt");
	leaveJoinFile << "#r/R join leave exchange" << endl;
	for(int i=0;i<rScaledJoin->GetXaxis()->GetNbins();i++)
	{
		leaveJoinFile << i*0.05 << " " << rScaledJoin->GetBinContent(i) << " " << rScaledLeave->GetBinContent(i) << " " << rExchange->GetBinContent(i) << endl;
	}
	leaveJoinFile.close();
	
	TFile* LeaveJoin = new TFile("leavejoin.root","RECREATE");
	rScaledJoin->Write("rScaledJoin");
	rScaledLeave->Write("rScaledLeave");
	rExchange->Write("rExchange");
	LeaveJoin->Close();
	
	delete cE,rScaledJoin,rScaledLeave,rExchange;
	delete hA,q;
	
	//Plot radii over time
	plotRadii(numSteps,bulkEdge,monoEdge,"radii");
// 	cout << "Bulk Mono:" << endl;
// 	for(int i=0;i<numFrames;i++)
// 		cout << bulkEdge[i] << " " << monoEdge[i] << endl;
// 	plotRadii(numSteps,bulkEdge,monoEdge,"a");
	cout << "Done!"<<endl;
	return 0;
}

TGraph *horizontalHist(TH1D* hist)
{
	char *title, *xLabel, *yLabel;
	int n;
	double *x,*y, xlo, xhi;

	//Copy points
	TGraph* g =	new TGraph(hist);
	n=g->GetN();
	x=g->GetX();
	y=g->GetY();
	xlo=hist->GetXaxis()->GetXmin();
	xhi=hist->GetXaxis()->GetXmax();

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
	g1->Draw("AL");
	g1->GetYaxis()->SetRangeUser(xlo,xhi);
	g1->Draw("AL");

	//Delete
	//delete g;
	return g1;
}

//Find the edge of droplet by tanh fitting for each row given TH2D
TGraph* findBoundaryPoints(TH2D* hist,char* aOrR,double& xBulkMax,double &xMonoEdge,int stepNum,int timestep)
{
	//Best guess for parameters
	double guess;
	
	//Whether the fit was Successful
	bool success;
	
	//Number of bins
	int nx=hist->GetNbinsX();
	int ny=hist->GetNbinsY();
	
	//Parameters: liquid density,interface width,center
	double ld,w,c;
	
	//Check for valid results
	double tmp;
	
	//Lists of boundary point coordinates
	double xCoords[nx+ny];
	double yCoords[nx+ny];
	
	//Number of rows containing the droplet
	int n=0;
	
	//Fitting tanh function
	TF1* tanhFit = new TF1("tanhFit","[0]/2*(1-tanh((x-[2])/(2*[1])))",0,200);
	
	//For each row
	for(int j=1;j<=ny;j++)
	{
		//Set Parameters based on counting guess
		guess=guessRowBoundary(hist,j);
		tanhFit->SetParameters(1,10,guess);
// 		cout<<"Row "<<j<<endl;
		//Create projection
		TH1D* px = hist->ProjectionX("px",j,j);
		
		//Create graph
		TGraph* gR = new TGraph(px);
		gR->Fit(tanhFit,"Q");

		//Get Parameters
		ld=tanhFit->GetParameter(0);
		w=tanhFit->GetParameter(1);
		c=tanhFit->GetParameter(2);
		
// 		//Save
// 		if(false)
// 		{
// 			stringstream cName,fName;
// 			cName << aOrR << "Row" << "Step" << setw(3) << setfill('0') <<  stepNum << "_" << setw(3) << setfill('0') << j;
// 			fName << "img/tanh/" << cName.str() << ".png";
// 			TCanvas* cR = new TCanvas();
// 			cR->cd();
// 			gR->SetTitle(cName.str().data());
// 			0SetMarkerStyle(20);
// 			gR->SetMinimum(0.0);
// 			gR->SetMaximum(1.5);
// 			gR->Draw("APL");
// 			cR->SaveAs(fName.str().data());
// 			delete cR;
// 		}
		
		delete px;
		
		//ld>0.5=>fit is valid and intersects 0.5=>row contains droplet
		if(ld>0.5)
		{
// 			cout << "Found" << endl;
			//Solve for x coordinate where the tanhFit(x)=0.5
			tmp=2*w*atanh(1-1/ld)+c;
			if(tmp>0)
			{
				xCoords[n]=tmp;
				yCoords[n]=hist->GetYaxis()->GetBinCenter(j);
			
				//Found another row containing droplet
				n++;
			}
			
		}

		delete gR;
	}
	
	//Choose the value from the second row to be the x coordinate after which to discard points for circle fitting to determine the bulk-monolayer interface
	xBulkMax=xCoords[1];

	//Choose the higher of the values from the first two rows to be the edge of the monolayer
// 	xMonoEdge=max(xCoords[0],xCoords[1]);
	xMonoEdge=xCoords[0]; //Use first row
	
	
	tanhFit->SetParameters(1,10,50);
	
	//For each column
	for(int i=1;i<=nx;i++)
	{
		//Best guess based on counting
		guess=guessColBoundary(hist,i);
		tanhFit->SetParameters(1,10,guess);
		
// 		cout << "Col "<<i<<endl;
		//Create projection
		TH1D* py = hist->ProjectionY("py",i,i);
		
		//Create graph
		TGraph* gC = new TGraph(py);
		gC->Fit(tanhFit,"Q");
		
		//Get Parameters
		ld=tanhFit->GetParameter(0);
		w=tanhFit->GetParameter(1);
		c=tanhFit->GetParameter(2);
		
		//Assume the fit is okay until proven guilty
		success=true;
		if(abs(c)>1000) success=false;
		
		delete py;
		
// 		Save
		
//		if(timestep==6048000)
//		{
//			stringstream cName,fName;
//			cName << aOrR << "Col" << "Step" << setw(3) << setfill('0') << timestep << "_" << setw(3) << setfill('0') << i;
//			system("mkdir -p img/tanh");
//			fName << "img/tanh/" << cName.str() << ".png";
//			TCanvas* cC = new TCanvas();
//			cC->cd();
//			gC->SetTitle(cName.str().data());
//			gC->SetMarkerStyle(20);
//			gC->SetMinimum(0.0);
//			gC->SetMaximum(1.5);
//			gC->Draw("APL");
//			cC->SaveAs(fName.str().data());
//			delete cC;
//		}
		
			
		//ld>0.5=>fit is valid and intersects 0.5=>column contains droplet
		//success => fit didn't blow up
		if(success && ld>0.5)
		{
// 			cout << "Found" << endl;
			//Solve for x coordinate where the tanhFit(x)=0.5
			tmp=2*w*atanh(1-1/ld)+c;
			if(tmp>=0)
			{
				xCoords[n]=hist->GetXaxis()->GetBinCenter(i);
				yCoords[n]=tmp;
				//Found another column containing droplet
				n++;
			}
		}

		delete gC;
	}
	
	delete tanhFit;
	
	//Overall 1D Graph of boundary points
	TGraph* pointsGraph = new TGraph(n,xCoords,yCoords);
	pointsGraph->SetName(aOrR);
	pointsGraph->SetTitle(aOrR);
	
// 	TCanvas* c5 = new TCanvas();
// 	c5->cd();
	pointsGraph->Draw("same PL");
// 	c5->SaveAs("test.png");
// 	delete c5;
	
	return pointsGraph;
}


//Given a graph of points, a maximum x value, and the y coordinate of the interface with the substrate, fit a circle to the points and find the intersection of the circle with the substrate interface. The result is the bulk-monolayer interface
CircleFit fitCircle(TGraph* g,double xMax,int stepNum,int timestep)
{
	//Fit circle and intersect it with a constant c (Interface)
	char* name1 = (char*) g->GetTitle();
	stringstream nameStream;
	nameStream << name1 << setw(8) << setfill('0') << timestep;
	char* name2 = (char*) nameStream.str().data();
	int n=g->GetN();
	vector<double> x(n),y(n);
	double xTest,yTest;
	
// 	//Limits for circle fitting for first 50 timesteps
// 	double lowLim=30;
// 	double highLim=50;
// 	
	//Number of points with x<xMax
	int nValid=0;
	
	for(int i=0;i<n;i++)
	{
		//Test the point before including it
		g->GetPoint(i,xTest,yTest);
		
		//Only use point if it is less than xMax
		if(xTest<=xMax)
		{
// 			//For the first 50 timesteps, only use about the bottom half of the droplet
// 			if(timestep%2000<50)
// 			{
// 				g->GetPoint(i,x[nValid],y[nValid]);
// 				nValid++;
// 			}
// 			
// 			else if( (lowLim<=yTest) && (yTest<highLim) )
// 			{
				g->GetPoint(i,x[nValid],y[nValid]);
				nValid++;
// 			}
		}
// 		cout << "("<<x[i]<<","<<y[i]<<")"<<endl;
	}
	x.resize(nValid);
	y.resize(nValid);
	
	CircleFit C(name2,x,y);
// 	system("mkdir -p img/circles");
	C.Draw();
	C.Print();
	
	return C;
	
}

//Plot bulk and mono radii over time
void plotRadii(int numSteps,double* bulkEdge,double* monoEdge,char* name)
{
	//Timestep array
	double t[numSteps];
	for(int i=0;i<numSteps;i++)
		t[i]=i;
	
// 	TFile* file = new TFile("radii.root","new");
	
	TCanvas* cRadii = new TCanvas();
	cRadii->cd();
	
	TGraph* gBulk = new TGraph(numSteps,t,bulkEdge);
	TGraph* gMono = new TGraph(numSteps,t,bulkEdge);
	
// 	gBulk->Write("gBulk");
// 	gMono->Write("gMono");
	
	TMultiGraph* gMulti = new TMultiGraph();
	gMulti->Add(gBulk,"Bulk");
	gMulti->Add(gMono,"Mono");
	
// 	gMulti->Write("gMulti");
	
	gMulti->SetTitle("Radii");
	gMulti->Draw("APL");
	
// 	cRadii->Write("cRadii");
	
	stringstream fname;
	fname << name << ".png";
	
	cRadii->SaveAs(fname.str().data());
}

//Count the number of water atoms in the first timestep
int countAtoms(ifstream &inFile)
{
    //Count number of atoms
    bool countFlag=true;
    string line;
    int numAtoms=0;
    
    //Ignore the first 3 lines
	for(int i=0;i<3;i++) inFile.ignore(256,'\n');
    
    while(countFlag)
    {
	getline(inFile,line);
	
	//Count until reaching a line containing "TIMESTEP"
	if(line.find("TIMESTEP")!=string::npos||inFile.eof())
	{
	    countFlag=false;
		numAtoms-=1; //Account for the blank line
	}
	else
	    numAtoms++;
    }
	
    //Unset eof flag (if set) & return to beginning of file
    inFile.clear();
    inFile.seekg(0,ios::beg);
    
	cout << "Counted " << numAtoms << " atoms." << endl;
	
    return numAtoms;
}

//Split a string into a string array of words
vector<string> strSplit(string str)
{
	int len=str.length();
	stringstream ss(str);
	int numWords=1;
	
	//Count number of words
	for(int ch=0;ch<len;ch++)
		if(isspace(str[ch]))
			numWords++;
	
	//Allocate array
	vector<string> arr(numWords);
	
	//Read string into array
	for(int i=0;i<len;i++)
	ss >> arr[i];
	
	return arr;
}

//Get index and coordinates from string array of words 
void strToData(double *coords,double *velocities,double &dipole,string line)
{
    string str;
	
	int index;
    //Split string into array of words
    vector<string> strArr=strSplit(line);
	
    //Get index
    str=strArr[0];
    index=atoi(str.data());
    
	//Save values
    //string -> cstring -> double
    for(int i=0;i<3;i++)
    {
		//Get coordinates from the second, third, and fourth elements in the array
		str=strArr[1+i];
		*coords++=atof(str.data());

		//Get velocities from the sixth, seventh, and eighth elements in the array
		str=strArr[5+i];
		*velocities++=atof(str.data());
    }
    
    //Save dipole moment cos(theta)
    dipole=atof(strArr[10].data());
}

//Find the mean of a double vector
double mean(vector<double> v)
{
	double m=0;
	for(vector<double>::iterator it = v.begin(); it != v.end(); ++it)
		m+=*it;
	m/=v.size();
	
	return m;
}

//Choose the highest value in a vector
double max(vector<double> v)
{
	int n=v.size();
	double max=v[0];
	
	for(int i=0;i<n;i++)
	{
		if(v[i]>max) 
			max=v[i];
	}
	
	return max;
		
}

//Square
inline double square(double x) {return x*x;}

//Arctanh
inline double atanh(double x) {return log((1+x)/(1-x))/2;}

//Guess boundary of water molecule by counting for a single row
double guessRowBoundary(TH2D* hist,int j)
{
	//Guess - static in case none is found for a particular row - use last guess
	static double guess=0;
	
	//Boundary is where density=0.5
	const double cutoff=0.5;
	
	//Number of bins
	int nx=hist->GetNbinsX();

	//Scan columns from right to left
	for(int i=nx;i>0;i--)
	{
		if(hist->GetBinContent(i,j)>=cutoff)
		{
			guess=hist->GetXaxis()->GetBinCenter(i);
			break;
		}
	}

	return guess;
}

//Guess boundary of water molecule by counting for a single column
double guessColBoundary(TH2D* hist,int i)
{
	//Guess - static in case none is found for a particular column - use last guess
	static double guess=0;
	
	//Boundary is where density=0.5
	const double cutoff=0.5;
	
	//Number of bins
	int ny=hist->GetNbinsY();

	//Scan rows from top to bottom
	for(int j=ny;j>0;j--)
	{
		if(hist->GetBinContent(i,j)>=cutoff)
			guess=hist->GetYaxis()->GetBinCenter(j);
	}

	return guess;
}

//Find boundary of water molecule for a single column
void guessColBoundary(TH2D* hist,double cutoff,vector<double> &xBounds,vector<double> &yBounds)
{
	//Was a boundary found for this column?
	bool found;

	//Number of bins
	int nx=hist->GetNbinsX();
	int ny=hist->GetNbinsY();

	//Scan columns from left to right
	for(int i=1;i<=nx;i++)
	{
		found=false;

		//Scan rows from top to bottom
		for(int j=ny;j>0;j--)
		{
			if(hist->GetBinContent(i,j)>=cutoff)
			{
				yBounds[i-1] = hist->GetYaxis()->GetBinCenter(j);
				
				found=true;
				break;
			}
		}

		//Use zero if this row has no boundary
		if(found==false)
		{
			yBounds[i-1]=0;
		/*
			if(i==1)
			{
				yBounds[i-1]=0;
			}
			else
			{
				yBounds[i-1]=yBounds[i-2];
			}
		*/
		}
		//X Bound is always predictable
		xBounds[i-1] = hist->GetXaxis()->GetBinCenter(i);
	}
}


//Calculate MSD - mean squared displacement
vector<double> calculateMSD(vector<double> xi,vector<double> yi,vector<double> zi,vector<double> x, vector<double> y,vector<double> z)
{
	vector<double> msd(4);
	
	int numAtoms=xi.size();
	
	double xSum=0;
	double ySum=0;
	double zSum=0;
	
	for(int i=0;i<numAtoms;i++)
	{
		xSum+=square(x[i]-xi[i])/numAtoms;
		ySum+=square(y[i]-yi[i])/numAtoms;
		zSum+=square(z[i]-zi[i])/numAtoms;
	}
	
	msd[0]=xSum;
	msd[1]=ySum;
	msd[2]=zSum;
	msd[3]=(xSum+ySum+zSum)/3;
	
	return msd;
}

//Count the number of atoms which have migrated from the bulk to the monolayer this timestep
double bulkMonoExchange(vector<double>z,double zInterface,int stepsPerFrame)
{
	//Average several timesteps
	static int stepCounter=0;
	stepCounter++;
	
	//Number of atoms in the simulation
	int numAtoms=z.size();
	
	//Number of atoms which have moved from the bulk to the monolayer this timestep
	static double exchange=0;
	//Final result to return
	double result;
	//Number of atoms which were previously in the monolayer
	static int prevCount=0;
	//Number of atoms which are currently in the monolayer
	int currCount=0;
	
	//Count atoms below the bulk-mono threshold
	for(int i=0;i<numAtoms;i++)
	{
		if(z[i]<zInterface)
			currCount++;
	}
	
	//The number which have moved to the monolayer is the difference between the number this timestep and last timestep divided by the time passed (ps)
	exchange=(currCount-prevCount)/(2*stepsPerFrame);
	
	//Save current count for next time
	prevCount=currCount;
	
	if(stepCounter%stepsPerFrame==0)
	{
		result=exchange;
		exchange=0;
		return result;
	}
}

//Keep track of which atoms join the monolayer, and save their radius scaled by the base radius to the vector rScaledJoin
void joinMonoPosition(vector<double> r,vector<double> z,double zInterface,double baseRadius,TH1D* rScaledJoin)
{
// 	cout << "Declare vars" << endl;
	//Number of atoms in the simulation
	int numAtoms=z.size();
	//List of atoms currently in the monolayer
	vector<int> monoList;
	//List of atoms in the monolayer last timestep
	static vector<int> prevMonoList;
	
// 	cout << "For loop" << endl;
	//Save ids of atoms below the bulk-mono threshold
	for(int i=0;i<numAtoms;i++)
	{
		if(z[i]<zInterface)
			monoList.push_back(i);
	}
// 	cout << "Size" << endl;
	//Number of atoms in the monolayer
	int numMonoAtoms=monoList.size();
// 	cout << "Check atoms" << endl;
	//For each atom in the monolayer, check whether it was there last step
	for(int i=0;i<numMonoAtoms;i++)
	{
		//If it's new to the monolayer, save it's scaled radial position
		if( !isIn(monoList[i],prevMonoList) )
			rScaledJoin->Fill(r[i]/baseRadius);
	}
	
	//Save monoList for next time
	prevMonoList=monoList;
}

//Keep track of which atoms leave the monolayer, and save their radius scaled by the base radius to the vector rScaledLeave
void leaveMonoPosition(vector<double> r,vector<double> z,double zInterface,double baseRadius,TH1D* rScaledLeave)
{
	//Number of atoms in the simulation
	int numAtoms=z.size();
	//List of atoms currently in the monolayer
	vector<int> monoList;
	monoList.clear();
	//List of atoms in the monolayer last timestep
	static vector<int> prevMonoList;
	
	//Save ids of atoms below the bulk-mono threshold
	for(int i=0;i<numAtoms;i++)
	{
		if(z[i]<zInterface)
			monoList.push_back(i);
	}
	
	//Number of atoms in the monolayer
	int numPrevMonoAtoms=prevMonoList.size();
	
	//For each atom in the monolayer last time, check whether it's still there
	for(int i=0;i<numPrevMonoAtoms;i++)
	{
		//If has just left the monolayer, save it's scaled radial position
		if( !isIn(prevMonoList[i],monoList) )
			rScaledLeave->Fill(r[i]/baseRadius);
	}
	
	//Save monoList for next time
	prevMonoList=monoList;
}

bool isLess(int a,int b) { return (a<b); }

//Check whether x is in v
bool isIn(int x, vector<int>v)
{
	return binary_search(v.begin(),v.end(),x,isLess);
}

double findMinimum(vector<double> v)
{
	double min=v[0];
	for(int i=0;i<v.size();i++)
	{
		if(v[i]<min) 
		{
			min=v[i];
		}
	}
	return min;
}

double findMaximum(vector<double> v)
{
	double max=v[0];
	for(int i=0;i<v.size();i++)
	{
		if(v[i]>max) 
		{
			max=v[i];
		}
	}
	return max;
}

//Count the number of timesteps
int countSteps(ifstream &inFile)
{
    string line;
    int numSteps=0;
    int lineNum=0;
    
    //Count number of timesteps
    while(getline(inFile,line))
	{
		if(line.find("TIMESTEP")!=string::npos)
			numSteps++;
	}
   
    //Unset eof flag (if set) & return to beginning of file
    inFile.clear();
    inFile.seekg(0,ios::beg);
    
	cout << "Counted " << numSteps << " timesteps." << endl;
    return numSteps;
}

//Find the lowest atom ID for a water molecule in this simulation
int lowestID(ifstream &inFile)
{
	int id;
	
	//Ignore the first 3 lines
	for(int i=0;i<3;i++)
		inFile.ignore(256,'\n');
	
	//Get ID
	inFile >> id; 
	
	//Return to beginning
	inFile.clear();
    inFile.seekg(0,ios::beg);
	
	return id;
}

//Test whether file exists
bool file_exists(const string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}
