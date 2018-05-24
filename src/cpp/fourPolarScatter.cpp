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
int main(int argc,char* argv[])
{
	//Command line arguments: inLoc,outLoc
	cout << "Received " << argc+1 << " command line arguments." << endl;

	char* inLoc=argv[1];
	char* outLoc=argv[2];

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
	
	//Histogram limits
	double rlo=0;
	double rhi=175;
	int nr[4]={25,50,75,200};
	
	double zlo=25;
	double zhi=100;
	int nz[4]={25,50,75,200};

	//Radius
	double r,lowBin,upBin;

	//Bin variables
	double binArea;
	double binHeight;
	double binVolume;

	//Skip first 2 lines
	for(int i=0;i<2;i++) inFile.ignore(256,'\n');
	lineNum+=2;
	
	//Create Canvas
	TCanvas *c1 = new TCanvas("c1","c1",1440,900);
	c1->Divide(2,2,0,0);
	//TCanvas *c2 = new TCanvas();

	//Create Histograms
	TH2D *h1 = new TH2D("h1","h1",nr[0],rlo,rhi,nz[0],zlo,zhi);
	TH2D *h2 = new TH2D("h2","h2",nr[1],rlo,rhi,nz[1],zlo,zhi);
	TH2D *h3 = new TH2D("h3","h3",nr[2],rlo,rhi,nz[2],zlo,zhi);
	TH2D *h4 = new TH2D("h4","h4",nr[3],rlo,rhi,nz[3],zlo,zhi);
	//TH2D *h2 = new TH2D("h2","\rho(x,y)",50,-300,300,50,-300,300);
	//h1->SetStats(0);
	double colzMin=0.0;
	double colzMax=0.05;

	h1->SetMinimum(colzMin);
	h1->SetMaximum(colzMax);
	h2->SetMinimum(colzMin);
	h2->SetMaximum(colzMax);
	h3->SetMinimum(colzMin);
	h3->SetMaximum(colzMax);
	h4->SetMinimum(colzMin);
	h4->SetMaximum(colzMax);

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

			lowBin=lowerBin(r,rlo,rhi,nr[0]);
			upBin=upperBin(r,rlo,rhi,nr[0]);
			binArea=PI*(square(upBin)-square(lowBin));
			binHeight=upperBin(z[i],zlo,zhi,nz[0])-lowerBin(z[i],zlo,zhi,nz[0]);
			binVolume=binArea*binHeight;
			h1->Fill(r,z[i],1.0/binVolume);

			lowBin=lowerBin(r,rlo,rhi,nr[1]);
			upBin=upperBin(r,rlo,rhi,nr[1]);
			binArea=PI*(square(upBin)-square(lowBin));
			binHeight=upperBin(z[i],zlo,zhi,nz[1])-lowerBin(z[i],zlo,zhi,nz[1]);
			binVolume=binArea*binHeight;
			h2->Fill(r,z[i],1.0/binVolume);

			lowBin=lowerBin(r,rlo,rhi,nr[2]);
			upBin=upperBin(r,rlo,rhi,nr[2]);
			binArea=PI*(square(upBin)-square(lowBin));
			binHeight=upperBin(z[i],zlo,zhi,nz[2])-lowerBin(z[i],zlo,zhi,nz[2]);
			binVolume=binArea*binHeight;
			h3->Fill(r,z[i],1.0/binVolume);

			lowBin=lowerBin(r,rlo,rhi,nr[3]);
			upBin=upperBin(r,rlo,rhi,nr[3]);
			binArea=PI*(square(upBin)-square(lowBin));
			binHeight=upperBin(z[i],zlo,zhi,nz[3])-lowerBin(z[i],zlo,zhi,nz[3]);
			binVolume=binArea*binHeight;
			h4->Fill(r,z[i],1.0/binVolume);

			//h2->Fill(x[i],y[i]);

		}
		//Title
		stringstream title;
		title << "Timestep " << timestep;
		//h1->SetTitle(title.str().data());
		//h1->SetTitle("z(r)");
		//h2->SetTitle("Real Location");

		//Draw
		c1->cd(1);
//		gPad->SetTickx(2);
		h1->Draw("colz");

		c1->cd(2);
//		gPad->SetTickx(2);
//		gPad->SetTicky(2);
		h2->GetYaxis()->SetLabelOffset(0.01);
		h2->Draw("colz");

		c1->cd(3);
		h3->Draw("colz");

		c1->cd(4);
//		gPad->SetTicky(2);
		h4->Draw("colz");


		//c1->cd();
		//h1->Draw("colz");
		//c2->cd();
		//h2->Draw("colz");
		//Save
		stringstream filename;
		filename << outLoc << "/img/step" << setw(7) << setfill('0') << timestep << ".jpg";
		stringstream filename2;
		//filename2 << outLoc << "/img/absloc.jpg";
		c1->SaveAs(filename.str().data());
		//c2->SaveAs(filename2.str().data());
		h1->Reset();
		h2->Reset();
		h3->Reset();
		h4->Reset();
		
		//Only first timestep
 		//loopFlag1=false;
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

