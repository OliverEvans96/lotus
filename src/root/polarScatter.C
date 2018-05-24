//Read each timestep and plot z(r)

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include "TCanvas.h"
#include "TH2D.h"

using namespace std;

const double PI = 3.141592653589793;

//Count the number of water atoms in the first timestep
int countAtoms(ifstream &inFile);

//Split a string into a string array of words
vector<string> strSplit(string str);

//Get index and coordinates from string array of words
void strToData(double *coords,string line);

//Find the mean of a double vector
double mean(vector<double> v);

//Square
double square(double x);

//Round to upper bin boundary
double upperBin(double x,double xlo,double xhi,int n);

//Round to lower bin boundary
double lowerBin(double x,double xlo,double xhi,int n);

//Read each timestep and plot 3d scatter plot
int polarScatter(char *inLoc,char *outLoc)
{
	//Command line arguments: inLoc,outLoc
	cout << "Received " << 2  << " command line arguments." << endl;

//	char* inLoc=argv[1];
//	char* outLoc=argv[2];

	cout << "Open Stream" << endl;
	//Input Stream
	ifstream inFile;
	cout << "Test" << inLoc << " :) " << inLoc << endl;
	inFile.open(inLoc);
	cout << "Test1" << endl;

	if (inFile.good())
		cout << "Successfully opened " << inLoc << endl;
	else
		cout << "Failed to open " << inLoc << endl;

	cout << "Variables" << endl;
	
	//Variables
	string line;
	string input;
	int atomNum=0;
	int lineNum=1;
	double coords[3];
	int timestep;
	bool loopFlag1=true;
	bool loopFlag2=true;
	int numAtoms=countAtoms(inFile);

	//Coordinates
	vector<double> x(numAtoms);
	vector<double> y(numAtoms);
	vector<double> z(numAtoms);
	
	//Bin size (A)
	double binSize=3.0;

	//Histogram limits
	double rlo=0;
	double rhi=175;
	int nr=(int)round((rhi-rlo)/binSize);
	
	double zlo=25;
	double zhi=100;
	int nz=(int)round((zhi-zlo)/binSize);

	//Radius
	double r,lowBin,upBin;

	//Bin variables
	double binArea;
	double binHeight;
	double binVolume;

	//Conversion factors
	double gramsPerWaterMolecule=2.9914e-23;
	double A3perCC=1.0e24;
	double convFact=gramsPerWaterMolecule*A3perCC;


	//Skip first 2 lines
	for(int i=0;i<2;i++) inFile.ignore(256,'\n');
	lineNum+=2;
	
	//Create Canvas
	TCanvas *c1 = new TCanvas("c1","c1",1440,900);

	//Create Histograms
	TH2D *h1 = new TH2D("h1","h1",nr,rlo,rhi,nz,zlo,zhi);
	h1->SetStats(0);

	double colzMin=0.0;
	double colzMax=1.0;

	h1->SetMinimum(colzMin);
	h1->SetMaximum(colzMax);

	//Read data
	while(!inFile.eof()&&loopFlag1)
	{
		
		//Read the header
		inFile >> line >> timestep;
		inFile.ignore();
		lineNum++;
		
		cout << "Timestep " << timestep << " @ line " << lineNum-1 << endl;
		
		loopFlag2=true;
		
		atomNum=0;

		//Read atoms
		while(loopFlag2)
		{
			//cout << "Test0" << endl;
			getline(inFile,line);
			//cout << "Test0.5" << endl;
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
		static double x0=mean(x);
		static double y0=mean(y);
		//cout << "Centroid: (" << x0 << "," << y0 << ")" << endl;
		
		//Fill Histograms
		for(int i=0;i<numAtoms;i++)
		{
			r=sqrt(square(x[i]-x0)+square(y[i]-y0));

			lowBin=lowerBin(r,rlo,rhi,nr);
			upBin=upperBin(r,rlo,rhi,nr);
			binArea=PI*(square(upBin)-square(lowBin));
			binHeight=upperBin(z[i],zlo,zhi,nz)-lowerBin(z[i],zlo,zhi,nz);
			binVolume=binArea*binHeight;
			h1->Fill(r,z[i],convFact/binVolume);

		}
		//Title
		stringstream title;
		title << "50A/atom1 Density (g/cc):  " << timestep;
		h1->SetTitle(title.str().data());

		//Axis labels
		h1->GetXaxis()->SetTitle("#rho (#AA)");
		h1->GetXaxis()->CenterTitle();
		h1->GetYaxis()->SetTitle("z (#AA)");
		h1->GetYaxis()->CenterTitle();

		//Draw
		c1->cd(1);
		h1->Draw("colz");

		//Save
		stringstream filename;
		filename << outLoc << "/img/step" << setw(7) << setfill('0') << timestep << ".jpg";
		stringstream filename2;
		c1->SaveAs(filename.str().data());
//		h1->Reset();
		
		//Only first timestep
 		loopFlag1=false;
	}
	
	inFile.close();

	return 0;
}

//Count the number of water atoms in the first timestep
int countAtoms(ifstream &inFile)
{
    //Count number of atoms
    bool countFlag=true;
    string line;
    int numAtoms=0;
	int numLines;	
    
    //Ignore the first 3 lines
	for(int i=0;i<3;i++) inFile.ignore(256,'\n');
    
    while(countFlag)
    {
	getline(inFile,line);
	
	//Count until reaching a line containing "TIMESTEP"
	if(line.find("TIMESTEP")!=string::npos||inFile.eof())
	{
	    countFlag=false;
		numLines-=1; //Account for the blank line
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

//Square
double square(double x) {return x*x;}

//Round to lower bin boundary
double lowerBin(double x,double xlo,double xhi,int n)
{
	return floor(n*(x-xlo)/(xhi-xlo))*(xhi-xlo)/n+xlo;
}

//Round to upper bin boundary
double upperBin(double x,double xlo,double xhi,int n)
{
	return ceil(n*(x-xlo)/(xhi-xlo))*(xhi-xlo)/n+xlo;
}

