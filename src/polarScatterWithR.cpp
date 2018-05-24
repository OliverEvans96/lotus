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
#include "TCanvas.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include "TFile.h"
#include "CircleFitClass.h"

using namespace std;

//There is an issue with the tanh fitting for some columns in the atom1


const double PI = 3.141592653589793;

//Count the number of water atoms in the first timestep
int countAtoms(ifstream &inFile);

//Split a string into a string array of words
vector<string> strSplit(string str);

//Get index and coordinates from string array of words
void strToData(double *coords,string line);

//Choose the higher of two doubles
double max(double a,double b);

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

//Find boundary points by tanh fitting for each row
TGraph* findBoundaryPoints(TH2D* hist,char* aOrR,double& xBulkMax,double &xMonoEdge,int stepNum,int timestep);

//Given a graph of points, a maximum x value, and the y coordinate of the interface with the substrate, fit a circle to the points and find the intersection of the circle with the substrate interface. The result is the bulk-monolayer interface
double findInterface(TGraph* g,double xMax,double c,int stepNum,int timestep);

//Plot bulk and mono radii over time
void plotRadii(int numSteps,double* bulkEdge,double* monoEdge,char* name);

//Find the minimum value of a vector
double findMinimum(vector<double> v);

//Count the number of timesteps
int countSteps(ifstream &inFile);

//Check whether a file exists
bool file_exists(const std::string& name);

//Read each timestep and plot 3d scatter plot
int main(int argc,char* argv[])
{
	//Command line arguments: inLoc,outLoc
	cout << "Received " << argc+1 << " command line arguments." << endl;

	char* inLoc=argv[1];
	char* outLoc=argv[2];

	cout << "Open Stream" << endl;
	//Input Stream
	ifstream inFile(inLoc);

	if (inFile.good())
		cout << "Successfully opened " << inLoc << endl;
	else
		cout << "Failed to open " << inLoc << endl;

	cout << "Variables" << endl;
	
	//Output file stream
	FILE* aRadii = fopen("aRadii.txt","w");
// 	fprintf(aRadii,"#t bulk mono\n");
	fflush(aRadii);
	
	//Variables
	string line;
	string input;
	int atomNum=0;
	int lineNum=1;
	double coords[3];
	int timestep; //Real step number for simulation
	int stepNum=0; //Step number relative to this file starting with 0
	int frameNum=0;
	bool loopFlag1=true;
	bool loopFlag2=true;
	
	ifstream testFile(inLoc);
	int numSteps=countSteps(testFile);
	int numAtoms=countAtoms(testFile);
	testFile.close();
	
	//Number of timesteps to average together to make a frame
	int stepsPerFrame=5;
	
	//Number of frames (collections of timesteps)
	//If numSteps%numFramesPerStep!=0, then the last frame will have fewer than the rest
	int numFrames=(int) ceil(numSteps/stepsPerFrame);

	//Coordinates
	vector<double> x(numAtoms);
	vector<double> y(numAtoms);
	vector<double> z(numAtoms);
	vector<double> r(numAtoms);
	
	//Histogram limits
// 	int nr=100;
// 	double rlo=0;
// 	double rhi=250;
	
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
// 	double alo=rlo;
// 	double ahi=PI*rhi*rhi;
	int nA=(int) round((ahi-alo)/dA);
	
	double rlo=sqrt(alo/PI);
	double rhi=sqrt(ahi/PI);
	
	//Bin volume
	double dV=dA*dz;
	
	//Check that dA, alo, and ahi are compatible
	if( !(fmod((ahi-alo),dA)==0) )
	{
		cout << "(Ahi-Alo) is not divisible by dA" << endl;
		return 1;
	}
	
	cout << "nA=" << nA << endl;
	
	//Radius
// 	double lowBin,upBin;

	//Bin variables
// 	double binArea;
// 	double binHeight;
// 	double binVolume;

	//Conversion factors
	double gramsPerWaterMolecule=2.9914e-23;
	double A3perCC=1.0e24;
	double convFact=gramsPerWaterMolecule*A3perCC;


	//Skip first 2 lines
	for(int i=0;i<2;i++) inFile.ignore(256,'\n');
	lineNum+=2;
	
	//Create Canvas
// 	TCanvas *c1 = new TCanvas("c1","c1",1440,900);
	
	//Fix canvas height & scale width to preserve aspect ratio
	int cH=600;
	int cW=(int) round(cH*(rhi-rlo)/(zhi-zlo));
	TCanvas *c2 = new TCanvas("c2","c2",cW,cH);

	//Values to use for equal-area bins
	double aVals[nA];
	double rVals[nA];
	for(int i=0;i<=nA;i++)
	{
		aVals[i]=i*dA;
		rVals[i]=sqrt(aVals[i]/PI);
	}
		
	
	//Create Histograms
// 	TH2D *h1 = new TH2D("h1","h1",nr,rlo,rhi,nz,zlo,zhi);
	TH2D *hA = new TH2D("hA","hA",nA,rVals,nz,zlo,zhi);
// 	h1->SetStats(0);
	hA->SetStats(0);
	
	//Title
	stringstream title;
	title << "Density (g/cc):  " << timestep;
	hA->SetTitle(title.str().data());

	//Axis labels
	hA->GetXaxis()->SetTitle("r (#AA)");
	hA->GetXaxis()->CenterTitle();
	hA->GetYaxis()->SetTitle("z (#AA)");
	hA->GetYaxis()->CenterTitle();

	//Z boundary between monolayer and bulk
	double zInterface = hA->GetYaxis()->GetBinCenter(3);
	
	
	double colzMin=0.0;
	double colzMax=1.5;

// 	h1->SetMinimum(colzMin);
// 	h1->SetMaximum(colzMax);
	hA->SetMinimum(colzMin);
	hA->SetMaximum(colzMax);

// 	Edge r coordinates
	double bulkEdgeR[numFrames];
	double monoEdgeR[numFrames];
	double bulkEdgeA[numFrames];
	double monoEdgeA[numFrames];
	
	//Set all values to zero
	for(int i=0;i<numSteps;i++)
	{
		bulkEdgeR[i]=0;
		monoEdgeR[i]=0;
		bulkEdgeA[i]=0;
		monoEdgeA[i]=0;
	}
	
	//Read data
	//For each timestep
	while(!inFile.eof()&&loopFlag1)
	{
		cout << "Step " << stepNum << endl;
		//Read the header
		inFile >> line >> timestep;
		inFile.ignore();
		lineNum++;
		
// 		cout << "Timestep " << timestep << " @ line " << lineNum-1 << endl;
		
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
				strToData(coords,line);
				x[atomNum]=coords[0];
				y[atomNum]=coords[1];
				z[atomNum]=coords[2];
				atomNum++;
			}
		}
		
		//Center
		//If the center has already been determined for this simulation, use that.
		//Otherwise, write it.
		double x0=mean(x);
		double y0=mean(y);
		
		if(file_exists("centroid.txt"))
		{
			ifstream centroidIn("centroid.txt");
			centroidIn >> x0 >> y0;
			centroidIn.close();
		}
		
		else
		{
			ofstream centroidOut("centroid.txt");
			centroidOut << x0 << y0;
			centroidOut.close();
		}
		
			
		//Fill Histograms
		for(int i=0;i<numAtoms;i++)
		{
			r[i]=sqrt(square(x[i]-x0)+square(y[i]-y0));

// 			lowBin=lowerBin(r[i],rlo,rhi,nr);
// 			upBin=upperBin(r[i],rlo,rhi,nr);
// 			binArea=PI*(square(upBin)-square(lowBin));
// 			binVolume=binArea*dz;
// 			
// 			h1->Fill(r[i],z[i],convFact/binVolume);
			hA->Fill(r[i],z[i],convFact/(dV*stepsPerFrame));

		}
		
		//Polar Scatter Plot
		ofstream psOut("pS.txt");
		for(int i=0;i<numAtoms;i++)
		{
			psOut << r[i] << " " << z[i] << endl;
			psOut.flush();
		}
		psOut.close();
		
		if( (stepNum!=0) && ((stepNum+1)%stepsPerFrame==0) )
		{
			
			//Draw
	// 		double xBulkMaxR,xMonoEdgeR;
	// 		c1->cd();
	// 		cout << "b1" << endl;
	// 		h1->Draw("colz");
	// 		TGraph* b1 = findBoundaryPoints(h1,"r",xBulkMaxR,xMonoEdgeR,stepNum);
	// 		b1->SetMarkerStyle(20);
	// 		b1->Draw("same PL");
	// 		stringstream rName;
	// 		rName << "img/colz/rbin" << setw(3) << setfill('0') << stepNum << ".jpg";
	// 		c1->SaveAs(rName.str().data());
	// 		
			double xBulkMaxA,xMonoEdgeA;
			c2->cd();
			hA->Draw("colz");
			TGraph* b2 = findBoundaryPoints(hA,"a",xBulkMaxA,xMonoEdgeA,stepNum,timestep);
			b2->SetMarkerStyle(20);
			b2->Draw("same L");
			stringstream aName;
			system("mkdir -p img/hist"); //Create directory if it doesn't exist
			aName << "img/hist/abin" << setw(3) << setfill('0') << stepNum << ".jpg";

			//Find bulk-monolayer interface
	// 		cout << endl;
	// 		cout << "Radius Binning" << endl;
	// 		double xBulkEdgeR=findInterface(b1,xBulkMaxR,zInterface,stepNum);
	// 		cout << "Bulk edge: " << xBulkEdgeR << endl;
	// 		bulkEdgeR[stepNum]=xBulkEdgeR;
	// 		cout << "Mono edge: " << xMonoEdgeR << endl;
	// 		monoEdgeR[stepNum]=xMonoEdgeR;
	// 		
			cout << endl;
			cout << "Area Binning" << endl;
			double xBulkEdgeA=findInterface(b2,xBulkMaxA,zInterface,stepNum,timestep);
			cout << "Bulk edge: " << xBulkEdgeA << endl;
			bulkEdgeA[frameNum]=xBulkEdgeA;
			cout << "Mono edge: " << xMonoEdgeA << endl;
			monoEdgeA[frameNum]=xMonoEdgeA;
			c2->SaveAs(aName.str().data());
			
			cout << endl;
			
			
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
			//Minimum
			double min=findMinimum(z);
			cout << "Minimum z: " << findMinimum(z) << endl;


			//Save
	// 		stringstream filename;
	// 		filename << outLoc << "/img/step" << setw(7) << setfill('0') << timestep << ".jpg";
	// 		c1->SaveAs(filename.str().data());
	// 		h1->Reset();
			hA->Reset();
			
			
			//Save radii to file
			fprintf(aRadii,"%03d %10.6f %10.6f \n",timestep,xBulkEdgeA,xMonoEdgeA);
			fflush(aRadii);
			
			frameNum++;
			
		}
		//Only first timestep
		if(stepNum==19)
			loopFlag1=false;
		
		//Increment relative timestep counter
		stepNum++;
	}
	
	inFile.close();
	fclose(aRadii);
	
	//Plot radii over time
	plotRadii(numSteps,bulkEdgeR,monoEdgeR,"r");
// 	cout << "Bulk Mono:" << endl;
// 	for(int i=0;i<numFrames;i++)
// 		cout << bulkEdgeA[i] << " " << monoEdgeA[i] << endl;
// 	plotRadii(numSteps,bulkEdgeA,monoEdgeA,"a");
	
	return 0;
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
// 			fName << "img/tanh/" << cName.str() << ".jpg";
// 			TCanvas* cR = new TCanvas();
// 			cR->cd();
// 			gR->SetTitle(cName.str().data());
// 			gR->SetMarkerStyle(20);
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
		
		if(i<20)
		{
			stringstream cName,fName;
			cName << aOrR << "Col" << "Step" << setw(3) << setfill('0') << stepNum << "_" << setw(3) << setfill('0') << i;
			system("mkdir -p img/tanh");
			fName << "img/tanh/" << cName.str() << ".jpg";
			TCanvas* cC = new TCanvas();
			cC->cd();
			gC->SetTitle(cName.str().data());
			gC->SetMarkerStyle(20);
			gC->SetMinimum(0.0);
			gC->SetMaximum(1.5);
			gC->Draw("APL");
			cC->SaveAs(fName.str().data());
			delete cC;
		}
		
			
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
// 	c5->SaveAs("test.jpg");
// 	delete c5;
	
	return pointsGraph;
}


//Given a graph of points, a maximum x value, and the y coordinate of the interface with the substrate, fit a circle to the points and find the intersection of the circle with the substrate interface. The result is the bulk-monolayer interface
double findInterface(TGraph* g,double xMax,double c,int stepNum,int timestep)
{
	//Fit circle and intersect it with a constant c (Interface)
	char* name1 = (char*) g->GetTitle();
	stringstream nameStream;
	nameStream << name1 << setw(3) << setfill('0') << stepNum;
	char* name2 = (char*) nameStream.str().data();
	int n=g->GetN();
	vector<double> x(n),y(n);
	double xTest,yTest;
	
	//Limits for circle fitting for first 50 timesteps
	double lowLim=30;
	double highLim=50;
	
	//Number of points with x<xMax
	int nValid=0;
	
	for(int i=0;i<n;i++)
	{
		//Test the point before including it
		g->GetPoint(i,xTest,yTest);
		
		//Only use point if it is less than xMax
		if(xTest<=xMax)
		{
			//For the first 50 timesteps, only use about the bottom half of the droplet
			if(timestep%2000<50)
			{
				g->GetPoint(i,x[nValid],y[nValid]);
				nValid++;
			}
			
			else if( (lowLim<=yTest) && (yTest<highLim) )
			{
				g->GetPoint(i,x[nValid],y[nValid]);
				nValid++;
			}
		}
// 		cout << "("<<x[i]<<","<<y[i]<<")"<<endl;
	}
	x.resize(nValid);
	y.resize(nValid);
	
	CircleFit C(name2,x,y);
	system("mkdir -p img/circles");
	C.Draw();
	C.Print();
	
	cout << "Intersect with " << c << endl;
	return C.Intersect(c);
	
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
	fname << name << ".jpg";
	
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
void strToData(double *coords,string line)
{
    string str;
	
	int index;
    //Split string into array of words
    vector<string> strArr=strSplit(line);
	
    //Get index
    str=strArr[0];
    index=atoi(str.data());
    
    //Get coordinates from the second, third, and fourth elements in the array
    //string -> cstring -> double
    for(int i=0;i<3;i++)
    {
		str=strArr[1+i];
		*coords++=atof(str.data());
    }
}

//Find the mean of a double vector
double mean(vector<double> v)
{
	double m=0;
	for(std::vector<double>::iterator it = v.begin(); it != v.end(); ++it)
		m+=*it;
	m/=v.size();
	
	return m;
}

//Choose the higher of two doubles
double max(double a,double b)
{
	if(a>b) return a;
	else return b;
}

//Square
double square(double x) {return x*x;}

//Arctanh
double atanh(double x) {return log((1+x)/(1-x))/2;}

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
    
    return numSteps;
}

//Test whether file exists
bool file_exists(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}
