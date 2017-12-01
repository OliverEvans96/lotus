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
#include "TGraph2D.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TText.h"
#include "TPaveText.h"

#include "CircleFitClass.h"
#include "Quiver/Quiver.h"
// #include "../src/FieldViz/FieldViz.h"

using namespace std;

const double PI = 3.141592653589793;

//Count the number of water atoms in the first timestep
int countAtoms(ifstream &inFile);

//Split a string into a string vector of words
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

//Given a TH1D and a two bin numbers, draw a line between the points and solve for where y=yc (y_cutoff)
double solveLinear(TH1D *hist,int bin1,int bin2,double yc);

//Guess boundary of water molecule by counting for a single column
void guessTanhFit(TH1D* hist,double* fitBounds,double &width,double &boundary,double &xlo,double &xhi);

//Draw a TH1D horizontally, returning a TGraph which should probably be deleted
TGraph *horizontalHist(TH1D* hist);

//Find boundary points by tanh fitting for each row
void findBoundaryPoints(TH2D* hist,TGraph *circlePointsGraph,char* aOrR,double *monoLimits,TH1D *hMono,double& rBulkMax,double &monoEdge,double &bulkEdge,int frameStep,TLegend* tanhLegend,TLine** tanhLines,TText **tanhTexts,TPaveText *tanhTextBox);
//void findBoundaryPoints(TH2D* hist,TGraph2D *circlePointsGraph,char* aOrR,double *monoLimits,TH1D *hMono,double& rBulkMax,double &monoEdge,int frameStep);

//Given a graph of points, a maximum x value, and the y coordinate of the interface with the substrate, fit a circle to the points
double fitCircle(TGraph* circlePointsGraph,CircleFit &circle,double xMax,int timestep);
//double fitCircle(TGraph2D* circlePointsGraph,CircleFit circle,double xMax,int timestep);

double solveTanhFit(TH1D* hist,TF1* tanhFit,double* fitBounds,int startBin,int frameStep,double bulkEdge,char *fitType,double pos,TLegend* tanhLegend,TLine** tanhLines,TText **tanhTexts,TPaveText *tanhTextBox);

//Calculate MSD - mean squared displacement
vector<double> calculateMSD(vector<double> xi,vector<double> yi,vector<double> zi,vector<double> x, vector<double> y,vector<double> z);

//Count the number of atoms which have migrated from the bulk to the monolayer this timestep
//double bulkMonoExchange(vector<double>z,double *monoLimits,int stepsPerFrame);

//Keep track of which atoms join the monolayer, and save their radius scaled by the base radius
int monoFlux(vector<double> r,vector<double>z,double* monoLimits,double baseRadius,TH1D* rScaledJoin,TH1D* rScaledLeave,int &nMono);

//Is a<b?
bool isLess(int a,int b);

//Is x in v?
bool isIn(int x, vector<int>v);

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

//round up to nearest multiple of a number
int roundUp(int numToRound, int multiple);

//Find z interface between monolayer and bulk
void findMonoLimits(TH1D *hWaterDens,double *monoLimits);

//Conversion factor - number density to water mass density
const double convFact=18/.60221409;

//Read each timestep and plot 3d scatter plot
int main(int argc,char* argv[])
{
    //Command line arguments: inLoc,outLoc
    cout << "Received " << argc-1 << " command line arguments." << endl;

    char* inLoc=argv[1];
    char* outLoc=argv[2];

    //Whether to run whole analysis or only find monolayer
    bool onlyFindInterface;
    cout << "argc=" << argc << endl;
    if(argc<4)
        onlyFindInterface=false;
    else
        onlyFindInterface=(strcmp(argv[3],"mono")==0);

    /*
    cout << "argv: " << endl;
    for(int i=0;i<argc;i++)
        cout << argv[i] << endl;
    cout << endl;
    */

    //z limits of the monolayer
    double monoLimits[2]={14.5,15.5};
    double monoWidth=monoLimits[1]-monoLimits[0];

    //While reading files, we don't yet know where the monolayer will be defined for this frame, so we save an atom's coordinates if it's in a broad range, then eliminate extras later.
    //double potentialMonoLimits[2]={14.5,15.5};
    //int nMonoAtoms=0;
    //int nPotentialMonoAtoms=0;
    //vector<double> potentialMonoR;
    //vector<double> potentialMonoZ;

    //Open monoLimits file
    /*
    fstream interfaceFile;
    if(!onlyFindInterface)
    {
        interfaceFile.open("../mono_limits.txt",ios::in);
    }
    else
        interfaceFile.open("../mono_limits.txt",ios::out);


    if (interfaceFile.good())
        cout << "Successfully opened " << "mono_limits.txt" << endl;
    else
        cout << "Failed to open " << "mono_limits.txt" << endl;

    if(!onlyFindInterface)
    {
        interfaceFile >> monoLimits[0] >> monoLimits[1];
        cout << "monoLimits=" << monoLimits[0] << " " << monoLimits[1] << endl;
    }
    */

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

    //Substrate density input file
    string dStr(argv[1]);
    stringstream dSS;
    dSS << dStr.substr(0,dStr.find("calculated.txt")) << "substrate_density_halfA.txt";
    cout << "Opening " << dSS.str().data() << endl;
    ifstream densFile(dSS.str().data());

    //Center of mass file
    ifstream centerFile("center_of_mass.txt");

    //Check files
    if (inFile.good())
        cout << "Successfully opened " << inLoc << endl;
    else
    {
        cout << "Failed to open " << inLoc << endl;
        return 1;
    }

    if (densFile.good())
        cout << "Successfully opened " << "substrate_density_halfA.txt" << endl;
    else
    {
        cout << "Failed to open " << "substrate_density_halfA.txt" << endl;
        return 1;
    }

    if(centerFile.good())
        cout << "Successfully opened center_of_mass.txt" << endl;
    else
    {
        cout << "Failed to open center_of_mass.txt" << endl;
        return 1;
    }

    //Substrate density variables
    int histDensBin;
    string densLine;
    double density;
    //Read first line (header comment)
    getline(densFile,densLine);

    //Split Header comment into vector of strings
    vector<string> densHeader = strSplit(densLine);

    //Count number of fields & allocate zVals array
    int nDensBinsFile = densHeader.size()-1;
    double *zVals = new double[nDensBinsFile];

    cout << "nDensBinsFile=" << nDensBinsFile << endl;
    //Remove "z=" from the beginning of each label and convert to decimal
    for(int i=0;i<nDensBinsFile;i++)
    {
        zVals[i]=atof(densHeader[i+1].substr(2).data()); //substr(2) means char 2 (counting from 0 - actually 3rd character) until the end
    }
    double dLoFile=zVals[0];
    double dZDens=zVals[1]-zVals[0];
    double dHiFile=zVals[nDensBinsFile-1]+dZDens;

    //Density histogram limits - dz for file & hist should be same, everything else can be different
    double dLoHist=0;
    double dHiHist=60;
    int nDensBinsHist=(int) floor((dHiHist-dLoHist)/dZDens);
    cout << "nDensBinsHist=" << nDensBinsHist << endl;

    //Ignore blank line
    getline(densFile,densLine);

    //Output file streams
    FILE* avgStepData = fopen("avgStepData.txt","w");
     fprintf(avgStepData,"%8.8s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s\n","#time","bulkEdge","monoEdge","contactAngle","dropletHeight","avgRadius","avgDipole","MSDx","MSDy","MSDz","MSDavg","frameFlux","dMe_dt","nMono","chi2s");
    fflush(avgStepData);

    FILE* instStepData = fopen("instStepData.txt","w");
     fprintf(instStepData,"%8.8s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s %15.15s \n","#time","avgRadius","avgDipole","MSDx","MSDy","MSDz","MSDavg","stepFlux");
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
    int stepNum=1; //Step number relative to this file starting with 1
    int frameNum=0;
    int frameStep; //First real timestep # of frame - used for naming, plots, etc.
    bool frameNamed=false; //Whether this frame has been named yet

    bool loopFlag1=true;
    bool loopFlag2=true;
    bool skipToEnd=false; //If true, skip first 485 steps

    //Find out some things about the file
    cout << endl;
    cout << "Analyzing input file" << endl;
    ifstream testFile(inLoc);
    int numSteps=countSteps(testFile);
    int numAtoms=countAtoms(testFile);
    int firstID=lowestID(testFile);
    testFile.close();

    //Number of timesteps to average together to make a frame
    int stepsPerFrame=5;

    //Is the number of steps divisible by the number of steps per frame?
    bool extraSteps = numSteps%stepsPerFrame;
    bool divisible = (extraSteps==0);

    //End of the penultimate frame (or middle of last frame if not divisible)
    int penultimateFrame = numSteps-extraSteps;

    cout << "extraSteps: " << extraSteps << endl;
    cout << "divisible: " << divisible << endl;
    cout << "penultimateFrame: " << penultimateFrame << endl;

    //Number of frames (collections of timesteps)
    //If not divisible, then the last frame will have more than the rest
    int numFrames=(int) floor(numSteps/stepsPerFrame);

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
    double x0,y0,z0;
    bool centroidKnown;

    //Radius of cylinder to use to calculate density(z) of water
    //Use half of the radius of the droplet
    int simPos=0; //When submitted by analyze.sh, whole path is included in argv[2]
    if(string(argv[2]).length()>40)
        simPos=44;
    double rDensCyl=atoi(string(argv[2]).substr(simPos,2).data())/2;
    cout << "RDENSCYL=" << rDensCyl << endl;

    //Step variables
    double rBulkMax=rDensCyl*2;
    double bulkEdge=rDensCyl*2;
    double tmpBulkEdge; //throwaway
    double monoEdge;
    double lastMonoEdge=0;
    double dropletHeight,contactAngle;
    double avgRadius,avgDipole;
    vector<double> msd;
    //double exchange;
    double min;
    int nMono;
    int bin1,bin2;

    //Flux for this timestep
    int stepFlux = 0;
    //Average flux for frame
    int frameFlux = 0;

    //Hist of where atoms join the monolayer
    TH1D* rScaledJoin = new TH1D("rScaledJoin","rScaledJoin",30,0,1.5);

    //Hist of where atoms leave the monolayer
    TH1D* rScaledLeave = new TH1D("rScaledLeave","rScaledLeave",30,0,1.5);

    double dz=1.5;
    double zlo=10;
    double zhi=142; //142 for 60A, 121 for others
    int nz=(int) round((zhi-zlo)/dz);

    //Check that dz, zlo, and zhi are compatible
    if( !(fmod((zhi-zlo),dz)==0) )
    {
        cout << "(zhi-zlo) is not divisible by dz" << endl;
        return 1;
    }

    cout << "nz=" << nz << endl;

    const double dA=500;
    const double alo=0;
    const double ahi=125500;
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

    //Skip first 2 lines
    for(int i=0;i<2;i++) inFile.ignore(256,'\n');
    lineNum+=2;

    //Fix canvas height & scale width to preserve aspect ratio
    //Added values are to account for window decorations which ROOT for some reason considers
    //Height and width must be divisible by 16 to play nice with ffmpeg
    int cH=1200+28;
    int cW=(int) roundUp((int)floor(cH*(rhi-rlo)/(zhi-zlo)),16)+4;
    cout << "cW=" << cW%16 << endl;
    TCanvas *cD = new TCanvas("cD","cD",cW,cH); //Dipole
    TCanvas *cA = new TCanvas("cA","cA",cW,cH); //2D density plot
    TCanvas *cQ = new TCanvas("cQ","cQ",cW,cH); //Quiver
    TCanvas *cVr = new TCanvas("cVr","cVr",cW,cH); //Plot v_r(z)
    TCanvas *cDens = new TCanvas("cDens","cDens",cW,cH); //Plot density of substrate and water

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
    TH1D *hMono = new TH1D("hMono","hMono",nA,rVals); //Only monolayer atoms - for calculating monoEdge
    //TH1D *hMonoProj;
    TH1D *hD = new TH1D("hD","hD",20,0,1);
    TH1D *hVr = new TH1D("hVr","hVr",nz,zlo,zhi); //For tracking v_r(z)
    TH1D *hZ = new TH1D("hZ","hZ",nz,zlo,zhi); //For counting # of atoms in each z bin
    TH1D *hWaterDens = new TH1D("hWaterDens","hWaterDens",nDensBinsHist,dLoHist,dHiHist);
    TH1D *hSubstrDens = new TH1D("hSubstrDens","hSubstrDens",nDensBinsHist,dLoHist,dHiHist);

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

    hWaterDens->GetXaxis()->SetTitle("z (#AA)");
    hWaterDens->GetXaxis()->CenterTitle();
    hWaterDens->GetYaxis()->SetTitle("#rho (g/cm^{3})");
    hWaterDens->GetYaxis()->CenterTitle();

    //Circle fit graph
    TGraph *circlePointsGraph = new TGraph();
    //TGraph2D *circlePointsGraph = new TGraph2D();

    //////////////////////
    // HIST ANNOTATIONS //
    //////////////////////

    //Bulk & Mono radius vertical lines
    TLine *bulkEdgeLine = new TLine(0,zlo,0,zhi); //Bulk edge
    TLine *monoEdgeLine = new TLine(0,zlo,0,zhi); //Mono edge
    TLine *heightLine = new TLine(0,zlo,rhi,zhi); //Droplet height
    TLine *monoHiLine = new TLine(0,zlo,rhi,zhi); //top of monolayer
    TLine *monoLoLine= new TLine(0,zlo,rhi,zhi); //bottom of monolayer

    //Density plot
    TLine *monoHiLineDens = new TLine(zlo,0,zhi,6); //top of monolayer
    TLine *monoLoLineDens = new TLine(zlo,0,zhi,6); //bottom of monolayer

    //Line properties
    bulkEdgeLine->SetLineWidth(3);
    bulkEdgeLine->SetLineColor(kGreen);
    monoEdgeLine->SetLineWidth(3);
    monoEdgeLine->SetLineColor(kRed);
    heightLine->SetLineWidth(3);
    heightLine->SetLineColor(kOrange+3);
    monoHiLine->SetLineWidth(3);
    monoHiLine->SetLineColor(kBlue);
    monoHiLine->SetLineStyle(9);
    monoLoLine->SetLineWidth(3);
    monoLoLine->SetLineColor(kBlue);
    monoLoLine->SetLineStyle(2);

    //Density Plot
    monoHiLineDens->SetLineWidth(3);
    monoHiLineDens->SetLineColor(kRed);
    monoLoLineDens->SetLineWidth(3);
    monoLoLineDens->SetLineColor(kGreen);

    //2D Density Hist Legend
    TLegend *hALegend = new TLegend(.65,.65,.85,.85);
    hALegend->AddEntry(circlePointsGraph,"Droplet boundary","lp");
    hALegend->AddEntry(bulkEdgeLine,"Bulk radius","l");
    hALegend->AddEntry(monoEdgeLine,"Mono radius","l");
    hALegend->AddEntry(heightLine,"Droplet height","l");
    hALegend->AddEntry(monoHiLine,"Mono top","l");
    hALegend->AddEntry(monoLoLine,"Mono bottom","l");

    //Text to show data
    TPaveText *textBox = new TPaveText();
    textBox->SetX1NDC(.65);
    textBox->SetY1NDC(.5);
    textBox->SetX2NDC(.85);
    textBox->SetY2NDC(.625);

    TText *cAText=textBox->AddText("Contact angle"); //Contact angle
    TText *dHText=textBox->AddText("Droplet height"); //Droplet height
    TText *bEText=textBox->AddText("Bulk edge"); //Bulk edge
    TText *mEText=textBox->AddText("Mono edge"); //Mono edge

    textBox->SetShadowColor(0);
    textBox->SetTextSize(0.025);

    //////////////////////
    // TANH ANNOTATIONS //
    //////////////////////

    //Lines (There are 7)
    TLine *ldTanhLine = new TLine(rlo,zlo,rhi,zlo); //Liquid density (horizontal)
    TLine *solTanhLine = new TLine(0,zlo,0,zhi); //Final solution
    TLine *x0TanhLine = new TLine(0,zlo,0,zhi); //x0 final
    TLine *x0guessTanhLine = new TLine(0,zlo,0,zhi); //x0 guess
    TLine *lowGuessTanhLine = new TLine(0,zlo,0,zhi); //Lower edge of boundary
    TLine *hiGuessTanhLine = new TLine(0,zlo,0,zhi); //Upper edge of boundary
    TLine *lowBinTanhLine = new TLine(0,zlo,0,zhi); //Lowest bin used

    //Line properties
    ldTanhLine->SetLineColor(kYellow);
    ldTanhLine->SetLineWidth(3);
    solTanhLine->SetLineColor(kBlack);
    solTanhLine->SetLineWidth(3);
    x0TanhLine->SetLineColor(kBlue);
    x0TanhLine->SetLineWidth(3);
    x0guessTanhLine->SetLineColor(kCyan);
    x0guessTanhLine->SetLineWidth(3);
    lowGuessTanhLine->SetLineColor(kViolet);
    lowGuessTanhLine->SetLineWidth(3);
    hiGuessTanhLine->SetLineColor(kGreen);
    hiGuessTanhLine->SetLineWidth(3);
    lowBinTanhLine->SetLineColor(kOrange+7);
    lowBinTanhLine->SetLineWidth(3);

    //Tanh Legend
    TLegend *tanhLegend = new TLegend(.65,.65,.85,.85);
    tanhLegend->AddEntry(ldTanhLine,"ld","l");
    tanhLegend->AddEntry(solTanhLine,"sol","l");
    tanhLegend->AddEntry(x0TanhLine,"x0","l");
    tanhLegend->AddEntry(x0guessTanhLine,"x0guess","l");
    tanhLegend->AddEntry(lowGuessTanhLine,"lowGuess","l");
    tanhLegend->AddEntry(hiGuessTanhLine,"hiGuess","l");
    tanhLegend->AddEntry(lowBinTanhLine,"lowBin","l");

    //Text to show data
    TPaveText *tanhTextBox = new TPaveText();
    tanhTextBox->SetX1NDC(.65);
    tanhTextBox->SetY1NDC(.45);
    tanhTextBox->SetX2NDC(.85);
    tanhTextBox->SetY2NDC(.625);

    //Text lines (There are 5)
    TText *tanhTexts[5];
    tanhTexts[0]=tanhTextBox->AddText("pos: ");
    tanhTexts[1]=tanhTextBox->AddText("w_guess: ");
    tanhTexts[2]=tanhTextBox->AddText("w: ");
    tanhTexts[3]=tanhTextBox->AddText("x0_guess: ");
    tanhTexts[4]=tanhTextBox->AddText("x0: ");

    tanhTextBox->SetShadowColor(0);
    tanhTextBox->SetTextSize(0.025);

    //Create array of TLines to pass to solveTanhFit
    TLine *tanhLines[7] = {ldTanhLine,solTanhLine,x0TanhLine,x0guessTanhLine,lowGuessTanhLine,hiGuessTanhLine,lowBinTanhLine};

    /////////////////////////////
    // FieldViz Initialization //
    /////////////////////////////

//    double ff_eps = 5;
//    int lic_iter = 2;
//    double lic_expnt = .5;
//
//    FieldViz flowField(640,480);
//    flowField.SetColormap("/home/oge1/lammps/sapphire/analysis/src/FieldViz/phase_cmap.txt");
//    flowField.SetLICParams(lic_iter, lic_expnt);
//    flowField.Init();
//    flowField.SetBounds(rlo, rhi, zlo, zhi);
//    flowField.SetNAtoms(numAtoms);

    //////////
    // Misc //
    //////////

    //hA Zaxis limits
    double colzMin=0.0;
    double colzMax=1.5;

    //Circle for fitting
    //Use "cut" as second argument for 0.25<dens<0.75 cutoff
    CircleFit Circle(hA);
    double chi2s; //Fitting error

    //Test for hWaterDens
    int nWaterDens=0;
    int nWatersBelowMono=0;

    //Rate of change of monoEdge w.r.t. time.
    double dMe_dt;

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
    system("mkdir -p img/densityHalfA");
    system("mkdir -p img/vr");
    system("mkdir -p img/vp");
    system("mkdir -p img/all");
    system("mkdir -p img/circles");
    system("mkdir -p img/tanh");
    system("mkdir -p img/lic");
    system("mkdir -p fields");

    //Skip to end if desired (step 486,timestep 970000)
    if(skipToEnd)
    {
        cout << "Skipping to end!" << endl;
        while(lineNum<8374013)
        {
            getline(inFile,line);
            lineNum++;
        }
        stepNum=486;
    }

    ///////////////
    // TIME LOOP //
    ///////////////

    //For each timestep
    while( (stepNum<=numSteps) && (loopFlag1) && (!inFile.eof()))
    {
        //Read the header
        inFile >> line >> timestep;
        inFile.ignore();
        lineNum++;

        cout << endl;
        cout << "Step # " << stepNum << endl;
        cout << "Timestep " << timestep << " @ line " << lineNum-1 << endl;

        loopFlag2=true;

        atomNum=0;

        //Name frame
        if(!frameNamed)
        {
            frameStep=timestep;
            frameNamed=true;
        }

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
        if( (!initWritten) && (stepNum==1) )
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

        //Test - why is the last timestep repeated in instStepData.txt?
        cout << "Last line read: " << lineNum << endl;

        //Center of mass

        //Read center of mass file
        //Ignore timestep
        centerFile.ignore(256,' ');
        //Get center of mass
        centerFile >> x0 >> y0 >> z0;
        //Go to next line
        centerFile.ignore(256,'\n');

        //Use droplet xy center and substrate z center
        x0=mean(x);
        y0=mean(y);

        //Report findings
        cout << "{x0,y0,z0}={" << x0 << "," << y0 << "," << z0 << "}" << endl;

        /*
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
        */

        //Minimum
        min=findMinimum(z)-z0;

        //Track mono atoms
        if(trackMonoAtoms)
        {
            for(int i=0;i<numMonoIDs;i++)
            {
                monoR=sqrt(square(x[monoIDs[i]]-x0)+square(y[monoIDs[i]]-y0));
                monoZ=z[monoIDs[i]]-z0;
                //Fill graphs with radial positions of mono atoms
                monoGraph1->SetPoint(i,monoR,monoZ);
                monoGraph2->SetPoint(i,monoR,monoZ);
            }
        }


        /*
        cout << "convFact=" << convFact << endl;
        cout << "dA=" << dA << endl;
        cout << "dV=" << dV << endl;
        */

        //Fill Histograms
        for(int i=0;i<numAtoms;i++)
        {
            //Subtract center of mass from x,y,z coordinates
            x[i]-=x0;
            y[i]-=y0;
            z[i]-=z0;

            //Calculate radius
            r[i]=sqrt(square(x[i])+square(y[i]));
            p[i]=sqrt(square(r[i])+square(z[i]-zlo));

            vr[i]=(x[i]*vx[i]+y[i]*vy[i])/r[i]; // Chain rule

            hA->Fill(r[i],z[i],convFact/dV);
            hD->Fill(dipole,1/numAtoms);
            hVr->Fill(z[i],vr[i]);
            hZ->Fill(z[i]);
            q->Fill(r[i],z[i],vr[i],vz[i]);

            //Fill water density histogram
            //water=18amu, convert units, divide by volume
            if(r[i]<=rDensCyl)
            {
                hWaterDens->Fill(z[i],convFact/(dZDens*PI*rDensCyl*rDensCyl));
                nWaterDens++;
            }

            //Fill monolayer histogram
            if((monoLimits[0]<=z[i]) && (z[i]<=monoLimits[1]))
            {
                hMono->Fill(r[i],convFact/(dA*monoWidth));
                //nPotentialMonoAtoms++;
                //potentialMonoR.push_back(r[i]);
                //potentialMonoZ.push_back(z[i]);
            }
            else if(z[i]<monoLimits[0])
            {
                cout << "MOLECULE BELOW MONO THRESHOLD AT Z=" << z[i] << endl;
                nWatersBelowMono++;
            }
        }

        if(nWatersBelowMono>0)
            cout << "Found " << nWatersBelowMono << " waters below monolayer in step " << timestep << endl;
        nWatersBelowMono=0;

        //Fill Density histogram
        //Ignore timestep
        densFile >> densLine;
        for(int i=0;i<nDensBinsFile;i++)
        {
            densFile >> density;
            histDensBin = (int) floor((zVals[i]-dLoHist)*nDensBinsHist/(dHiHist-dLoHist))+1;
            //cout << "zVals[" << i << "]=" << zVals[i] << endl;
            //cout << (zVals[i]-dLoFile)*nDensBinsHist/(dHiHist-dLoHist) << endl;
            //cout << "Density " << zVals[i] << ": " << histDensBin << " = " << density << endl;
            hSubstrDens->SetBinContent(histDensBin,density);
        }
        //Go to next line
        densFile.ignore(256,'\n');

        //INST CALCULATIONS

        //Calculate average radius
        avgRadius=mean(r);

        //Calculate average dipole moment angle
        avgDipole=mean(cosT);

        //Calculate MSD
        msd=calculateMSD(xi,yi,zi,x,y,z);

        //Update lists of where atoms joined and left the monolayer
        stepFlux = monoFlux(r,z,monoLimits,bulkEdge,rScaledJoin,rScaledLeave,nMono);

        cout << "Step Flux: " << stepFlux << endl;

        //Calculate average flux for frame (gets reset after written every frame)
        frameFlux += stepFlux;

        //WRITE INST
        fprintf(instStepData,"%08d %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15d\n",timestep,avgRadius,avgDipole,msd[0],msd[1],msd[2],msd[3],stepFlux);

        ///////////////////////
        // FieldViz add data //
        ///////////////////////

//        flowField.SetData(&r[0], &z[0], &vr[0], &vz[0], numAtoms);
//        flowField.AddFieldFromData(ff_eps);

        //FRAME CALCULATIONS (5 timesteps or last step)

        //The last frame may contain more timesteps than other frames (generally 6, up to 9).
        //So if numSteps is not divisible by stepsPerFrame, then the following should not run for what would be the penultimate frame, and it should be combined with the remaining steps to form the final frame.
       //stepsPerFrame may need to be altered for the last timestep to correct some calculations.
        if( ((stepNum%stepsPerFrame==0)&&!(!divisible&&(stepNum==penultimateFrame))) || (stepNum==numSteps&&!divisible) )
        {
            cout << endl;
            cout << "///////////////" << endl;
            cout << "// FRAME " << setw(3) << stepNum/stepsPerFrame << " //" << endl;
            cout << "///////////////" << endl;
            cout << endl;
            cout << "frameStep: " << frameStep << endl;
            //stepsPerFrame needs to be altered on the last frame if it has extra steps
            if(!divisible && stepNum==numSteps)
                stepsPerFrame+=extraSteps;
            cout << stepsPerFrame << " steps this frame." << endl;
            cout << endl;

            //FRAME CALCULATIONS

            //Scale by stepsPerFrame since this is a frame averege
            hWaterDens->Scale(1.0/stepsPerFrame);
            hA->Scale(1.0/stepsPerFrame);
            hMono->Scale(1.0/stepsPerFrame);

            //Find monolayer
            //findMonoLimits(hWaterDens,monoLimits);

            /*
            //Write monolayer to file if appropriate
            if(onlyFindInterface)
            {
                interfaceFile << monoLimits[0] << " " << monoLimits[1];
                loopFlag1=false;
            }
            */

            /*
            //Determine which atoms are really in the monolayer
            for(int i=0;i<nPotentialMonoAtoms;i++)
            {
                if((monoLimits[0]<=potentialMonoZ[i]) && (z[i]<=potentialMonoR[1]))
                {
                    hMono->Fill(potentialMonoR[i],convFact/(dA*monoWidth));
                    nMonoAtoms++;
                }
            }
            */

            //Report findings
            //cout << "Found " << nPotentialMonoAtoms << " potentially in monolayer" << endl;
            //cout << " Of them, " << nMonoAtoms << " are really in the monolayer." << endl;

            /*
            //Reset potential monolayer variables
            nMonoAtoms=0;
            nPotentialMonoAtoms=0;
            potentialMonoR.clear();
            potentialMonoZ.clear();
            */

            cout << "Find Boundary" << endl;

            //Do Tanh fitting, find monolayer radius
            findBoundaryPoints(hA,circlePointsGraph,"a",monoLimits,hMono,rBulkMax,monoEdge,bulkEdge,frameStep,tanhLegend,tanhLines,tanhTexts,tanhTextBox);
            cout << "Mono edge: " << monoEdge << endl;

            cout << "Fit Circle" << endl;
            //Fit circle
            cA->cd();
            chi2s=fitCircle(circlePointsGraph,Circle,rBulkMax,frameStep);
            cout << "chi2s_ps:" << chi2s << endl;

            cout << "Find radii" << endl;

            //Find bulk radius
            // Previously using value from first bulk row
            // Now, just throw this away and use the tanhfit val from
            // first bulk row, which has already been found
            // in findBoundaryPoints
            tmpBulkEdge=Circle.Intersect(monoLimits[1]);
            cout << "Bulk edge: " << bulkEdge << endl;

            //Find contact angle
            cout << "Intersect with " << monoLimits[1] << endl;
            contactAngle=Circle.ContactAngle();

            //Find droplet height by intersecting circle with y-axis
            dropletHeight=Circle.Height();

            //Test number of molecules in hWaterDens
            cout << "STEP " << timestep << ": " << nWaterDens << " molecules in hWaterDens" << endl;
            nWaterDens=0;

            //PLOT TILES
            stringstream title;
            title << string(argv[2]).substr(simPos) << ": "  << frameStep;
            hA->SetTitle(title.str().data());

            title.str("");
            title << "Velocity Field: " << frameStep;
            q->SetTitle(title.str().data());

            title.str("");
            title << "Density: " << frameStep;
            hWaterDens->SetTitle(title.str().data());

            title.str("");
            title << "Dipole Moment Distribution: " << frameStep;
            hD->SetTitle(title.str().data());

            title.str("");
            title << "2D Radial Velocity: " << frameStep;
            hVr->SetTitle(title.str().data());

            //FRAME PLOTS

            //Draw dipole Histogram
            cD->cd();
            hD->Draw();
            title.str("");
            title << "img/dipole/step" << setw(8) << setfill('0') << frameStep << ".png";
            cD->SaveAs(title.str().data());

            //Quiver
            title.str("");
            title << "img/quiver/step" << setw(8) << setfill('0') << frameStep << ".png";
            q->Draw(cQ);
            cout << "Saving quiver" << endl;
            q->SaveAs(title.str().data());

            //Density plot
            title.str("");
            title << "img/densityHalfA/step" << setw(8) << setfill('0') << frameStep << ".png";
            cDens->cd();
            //Drawing options
            hWaterDens->SetLineColor(kBlue);
            hSubstrDens->SetLineColor(kOrange+3); //Brown
            hWaterDens->SetLineWidth(2);
            hSubstrDens->SetLineWidth(2);
            hWaterDens->Draw();
            hSubstrDens->Draw("SAME"); //Same canvas
            //Use OpenGL for antialiasing
            gStyle->SetCanvasPreferGL(true);
            //Axis limits
            hWaterDens->GetYaxis()->SetRangeUser(0,6);
            //Monolayer edges
            monoHiLineDens->SetX1(monoLimits[1]);
            monoHiLineDens->SetX2(monoLimits[1]);
            monoHiLineDens->Draw();
            monoLoLineDens->SetX1(monoLimits[0]);
            monoLoLineDens->SetX2(monoLimits[0]);
            monoLoLineDens->Draw();
            //Legend
            TLegend *densLeg = new TLegend(.75,.75,.85,.85);
            densLeg->AddEntry(hWaterDens,"Water");
            densLeg->AddEntry(hSubstrDens,"Substrate");
            densLeg->AddEntry(monoLoLineDens,"Mono lower limit","l");
            densLeg->AddEntry(monoHiLineDens,"Mono upper limit","l");
            densLeg->Draw();
            //Save
            cDens->SaveAs(title.str().data());

            //2D Velocity Plot
            //Scale values to create average
            for(int i=0;i<nz;i++)
            {
                sum=hVr->GetBinContent(i+1);
                num=hZ->GetBinContent(i+1);
                if(num!=0)
                    hVr->SetBinContent(i+1,abs(sum/num));
            //    cout << "hVr[" << i+1 << "] = " << hVr->GetBinContent(i+1) << endl;
            }
            //Plot
            cVr->cd();
            TGraph *g = horizontalHist(hVr);
            g->GetXaxis()->SetLimits(1e-6,1e-3);
            cVr->SetLogy();
            g->Draw("AP");
            cVr->Update();
            q->Draw(cVr);
            title.str("");
            title << "img/vr/step" << setw(8) << setfill('0') << frameStep << ".png";
            cVr->SaveAs(title.str().data());

            //Draw 2D density histogram
            cA->cd();
            hA->SetMinimum(colzMin);
            hA->SetMaximum(colzMax);
            hA->Draw("colz");
            //Draw circle
            TEllipse *circleEllipse = Circle.Draw();
            //Draw tangent line
            TLine *tangentLine = Circle.DrawTangentLine();
            //Add lines
            bulkEdgeLine->SetX1(bulkEdge);
            bulkEdgeLine->SetX2(bulkEdge);
            bulkEdgeLine->Draw();
            monoEdgeLine->SetX1(monoEdge);
            monoEdgeLine->SetX2(monoEdge);
            monoEdgeLine->Draw();
            heightLine->SetY1(dropletHeight);
            heightLine->SetY2(dropletHeight);
            heightLine->Draw();
            monoHiLine->SetY1(monoLimits[1]);
            monoHiLine->SetY2(monoLimits[1]);
            monoHiLine->Draw();
            monoLoLine->SetY1(monoLimits[0]);
            monoLoLine->SetY2(monoLimits[0]);
            monoLoLine->Draw();
            //Draw circle points graph
            circlePointsGraph->SetMarkerStyle(20);
            circlePointsGraph->Draw("same P");
            //Add legend (once)
            if(frameNum==0)
            {
                hALegend->AddEntry(tangentLine,"Tangent line","l");
            }
            hALegend->Draw();
            //Add data text box
            title.str("");
            title << "Contact angle: " << contactAngle;
            cAText->SetText(0,0,title.str().data());
            title.str("");
            title << "Droplet height: " << dropletHeight;
            dHText->SetText(0,0,title.str().data());
            title.str("");
            title << "Bulk radius: " << bulkEdge;
            bEText->SetText(0,0,title.str().data());
            title.str("");
            title << "Mono radius: " << monoEdge;
            mEText->SetText(0,0,title.str().data());
            textBox->Draw();

            //Name & Save
            stringstream aName;
            aName << "img/hist/abin" << setw(8) << setfill('0') << frameStep << ".png";
            cA->SaveAs(aName.str().data());

            //Delete ellipse and tangent line
            delete circleEllipse;
            delete tangentLine;

            //Draw mono atoms
            if(trackMonoAtoms)
            {
                monoGraph1->Draw("P");
                monoGraph2->Draw("P");

                stringstream mName;
                mName << "img/mono/monoHist" << setw(8) << setfill('0') << frameStep << ".png";
                cA->SaveAs(mName.str().data());
            }

            //Draw all together
            cAll->cd(1);
            gPad->Clear();
            cA->DrawClonePad();

            cAll->cd(2);
            gPad->Clear();
            cVr->DrawClonePad();

            cAll->cd(3);
            gPad->Clear();
            cQ->DrawClonePad();

            title.str("");
            title << "img/all/step" << setw(8) << setfill('0') << frameStep << ".png";

            //FOR SOME REASON, THE PROGRAM CRASHES HERE IF THE FILE IS SAVED
            //cAll->SaveAs(title.str().data());

            cout << "Reset histograms" << endl;
            //Reset histograms
            hA->Reset();
            hMono->Reset();
            hD->Reset();
            q->Reset();
            hWaterDens->Reset();
            hSubstrDens->Reset();

            //Calculate rate of change of monoEdge w.r.t. time (units A/step)
            dMe_dt=(monoEdge-lastMonoEdge)/stepsPerFrame;

            cout << endl;
            cout << "Frame Flux: " << frameFlux << endl;
            cout << "dMe_dt: " << dMe_dt << endl;
            cout << endl;

            //////////////////
            // FieldViz LIC //
            //////////////////

			/*
            title.str("");
            title << "img/lic/lic" << setw(8) << setfill('0') << frameStep << ".ppm";
			stringstream xFile, yFile;
			xFile << "fields/xfield" << setw(8) << setfill('0') << frameStep << ".txt";
			yFile << "fields/yfield" << setw(8) << setfill('0') << frameStep << ".txt";

            flowField.WriteField(xFile.str().data(),yFile.str().data(),"xy");
            flowField.PreProcess();
            flowField.PerformLIC();
            flowField.PostProcess();
            flowField.WriteImage2PPM(title.str().data());
            flowField.ResetField();
			*/

            //Save average data for frame to file
            fprintf(avgStepData,"%08d %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15d %15.6f %15d %15f\n",frameStep,bulkEdge,monoEdge,contactAngle,dropletHeight,avgRadius,avgDipole,msd[0],msd[1],msd[2],msd[3],frameFlux,dMe_dt,nMono,chi2s);
            fflush(avgStepData);
            cout << "Wrote step " << frameStep << " to avgStepData.txt" << endl;

            //Save monolayer radius
            lastMonoEdge=monoEdge;

            //Reset frame flux
            frameFlux = 0;
            frameNum++;

            //Reset name for next frame
            frameNamed = false;

            //Only first frame
            //loopFlag1=false;
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
    TH1D* rExchange = (TH1D*) rScaledJoin->Clone();
    //Normalize histograms
    rExchange->Add(rScaledLeave,-1);
    I=rScaledJoin->Integral();
    rScaledJoin->Scale(1/I);
    I=rScaledLeave->Integral();
    rScaledLeave->Scale(1/I);
    I=rExchange->Integral();
    rExchange->Scale(1/I);
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
        leaveJoinFile << i*(rScaledJoin->GetBinCenter(2)-rScaledJoin->GetBinCenter(1)) << " " << rScaledJoin->GetBinContent(i+1) << " " << rScaledLeave->GetBinContent(i+1) << " " << rExchange->GetBinContent(i+1) << endl;
    }
    leaveJoinFile.close();

    TFile* LeaveJoin = new TFile("leavejoin.root","RECREATE");
    rScaledJoin->Write("rScaledJoin");
    rScaledLeave->Write("rScaledLeave");
    rExchange->Write("rExchange");
    LeaveJoin->Close();

    delete cE,rScaledJoin,rScaledLeave,rExchange;
    delete hA,q,hD;
    delete hWaterDens,hSubstrDens;
    delete [] zVals;
    delete circlePointsGraph;

    cout << "Done!"<<endl;
    return 0;
}

TGraph *horizontalHist(TH1D* hist)
{
    char *title, *xLabel, *yLabel;
    int n;
    double *x,*y, xlo, xhi;

    //Copy points
    TGraph* g =    new TGraph(hist);
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
    delete g;
    return g1;
}

//Find the edge of droplet by tanh fitting for each row given TH2D
void findBoundaryPoints(TH2D* hist,TGraph *circlePointsGraph,char* aOrR,double *monoLimits,TH1D *hMono,double& rBulkMax,double &monoEdge,double &bulkEdge,int frameStep,TLegend* tanhLegend,TLine** tanhLines,TText **tanhTexts,TPaveText *tanhTextBox)
{
    cout << endl;
    cout << "FINDBOUNDARYPOINTS" << endl;

    //Best guess for parameters
    double guess;

    //Number of bins
    int nx=hist->GetNbinsX();
    int ny=hist->GetNbinsY();

    //Parameters: liquid density,interface width,center
    double ld,w,c;

    //Temporary values returned from solveTanhFit before saving
    double tmpx,tmpy;

    // bulkEdge, to be updated at the end of this function.
    double tmpBulkEdge;

    //Number of rows+columns containing the droplet
    int n=0;

    //Reset points graph
    circlePointsGraph->Set(0);
    //circlePointsGraph->Set(nx+ny);

    //Fitting tanh function
    TF1* tanhFit = new TF1("tanhFit","[0]/2*(1-tanh(4*(x-[2])/([1])))",0,300);

    //Set Bounds on parameters
    double fitBounds[6]={0.2,2.0,2,20,0,300};
    tanhFit->SetParLimits(0,fitBounds[0],fitBounds[1]); //ld
    tanhFit->SetParLimits(1,fitBounds[2],fitBounds[3]); //w
    tanhFit->SetParLimits(2,fitBounds[4],fitBounds[5]); //x0

    //Projections of hA
    TH1D *px,*py;

    //Bin center on axis perpendicular to projection
    double center=0;


    //Number of first z bin of hA above the monolayer
    py=hist->ProjectionY("py",1,1);
    int nz=py->GetNbinsX();
    double dz=py->GetBinWidth(1);
    double zlo=py->GetBinLowEdge(1);
    double zhi=py->GetBinLowEdge(nz)+dz;
    cout << "ZLO: " << zlo << " ZHI: " << zhi << endl;
    int firstBulkBin = (int) ceil((monoLimits[1]-zlo)*nz/(zhi-zlo))+1;
    cout << "FirstBulkBin=" << firstBulkBin << " : " << py->GetBinLowEdge(firstBulkBin) << endl;

    //Find monolayer edge
    //Fit hMono with tanh, find where density=0.5
    //Here we're passing monoEdge instead of bulkEdge because we're trying to find the monolayer, not bulk, radius.
    cout << "mono tanh fit" << endl;
    monoEdge=solveTanhFit(hMono,tanhFit,fitBounds,1,frameStep,monoEdge,"mono",center,tanhLegend,tanhLines,tanhTexts,tanhTextBox);

    //Only update monoEdge if point is good. Otherwise, use previous
    /*
    if(tmpx>0)
    {
        monoEdge=tmpx;
    }
    else
    {
        cout << "Mono Fit failed!!!" << endl;
        monoEdge=0;
    }*/
    cout << "Mono Radius: " << monoEdge << endl;

    //Row and column tanh fitting

    cout << "Row fits" << endl;
    //For each row
    for(int j=firstBulkBin+1;j<=ny;j++)
    {
        //Create projection
        px = hist->ProjectionX("px",j,j);
        center = hist->GetYaxis()->GetBinCenter(j);

        //Solve for x coordinate where tanhFit(x)=0.5
        //Include all bins
        tmpx=solveTanhFit(px,tanhFit,fitBounds,1,frameStep,bulkEdge,"row",center,tanhLegend,tanhLines,tanhTexts,tanhTextBox);
        //cout << "Row " << j << ": tmpx=" << tmpx << endl;
        if(tmpx>0)
        {
            //Found another row containing droplet
            tmpy=hist->GetYaxis()->GetBinCenter(j);
            circlePointsGraph->SetPoint(n,tmpx,tmpy);

            //Choose the value from the row above the monolayer to be the x coordinate after which to discard points for circle fitting to determine the bulk-monolayer interface
            if(j==firstBulkBin)
                rBulkMax=tmpx;
            n++;

            // If this is the first bulk row, save this value as bulkEdge
            if(j == firstBulkBin + 1)
                tmpBulkEdge = tmpx;

        }
    }

    cout << "column fits" << endl;
    //For each column
    for(int i=1;i<=nx;i++)
    {
        //Create projection
        py = hist->ProjectionY("py",i,i);
        center = hist->GetXaxis()->GetBinCenter(i);

        //Solve for y coordinate where the tanhFit(y)=0.5
        //Ignore first few bins in monolayer
        tmpy=solveTanhFit(py,tanhFit,fitBounds,firstBulkBin+1,frameStep,bulkEdge,"col",center,tanhLegend,tanhLines,tanhTexts,tanhTextBox);
        //cout << "Col " << i << ": tmpy=" << tmpy << endl;
        if(tmpy>=0)
        {
            //Found another column containing droplet
            tmpx=hist->GetXaxis()->GetBinCenter(i);
            circlePointsGraph->SetPoint(n,tmpx,tmpy);
            n++;
        }
    }
    cout << n << " points on graph" << endl;

    //Update graph properties
    circlePointsGraph->Set(n);
    circlePointsGraph->Sort();
    circlePointsGraph->SetName(aOrR);
    circlePointsGraph->SetTitle(aOrR);

    //cout << endl << "circlePointsGraph:" << endl;
    /*
    for(int i=0;i<n;i++)
    {
        //circlePointsGraph->GetPoint(i,tmpx,tmpy);
        tmpx=circlePointsGraph->GetX()[i];
        tmpy=circlePointsGraph->GetY()[i];
        //cout << "    (" << tmpx << "," << tmpy << ")" << endl;
    }
    */


    // Update bulkEdge
    bulkEdge = tmpBulkEdge;

    TCanvas* c5 = new TCanvas();
     c5->cd();
    circlePointsGraph->Draw("APL");
     c5->SaveAs("test.C");
     delete c5;

    delete tanhFit;
}


//Given a graph of points, a maximum x value, and the y coordinate of the interface with the substrate, fit a circle to the points and find the intersection of the circle with the substrate interface. The result is the bulk-monolayer interface
double fitCircle(TGraph* g,CircleFit &circle,double xMax,int timestep)
{
    //g->Draw("SAME");
    //Fit circle and intersect it with a constant c (Interface)
    char* name1 = (char*) g->GetTitle();
    stringstream nameStream;
    nameStream << name1 << setw(8) << setfill('0') << timestep;
    char* name2 = (char*) nameStream.str().data();
    int n=g->GetN();
    vector<double> x(n),y(n);
    double xTest,yTest;

    cout << "xMax=" << xMax << endl;

    //Sort points
    g->Sort();

//     //Limits for circle fitting for first 50 timesteps
//     double lowLim=30;
//     double highLim=50;
//
    //Number of points with x<xMax
    int nValid=0;

    for(int i=0;i<n;i++)
    {
        //Test the point before including it
        //g->GetPoint(i,xTest,yTest);
        xTest=g->GetX()[i];
        yTest=g->GetY()[i];

        //Only use point if it is less than xMax
        // I AM NOW USING ALL POINTS SINCE BAD ONES SHOULD ALREADY
        //if(xTest<=xMax)
        if(true)
        {
            x[nValid]=xTest;
            y[nValid]=yTest;
            //cout << x[nValid]<<" "<<y[nValid]<<endl;
            nValid++;
        }
        //Otherwise, mark the point as bad
        else
            circle.AddBadPoint(xTest,yTest);
    }
    x.resize(nValid);
    y.resize(nValid);

    //cout << "x-mean = " << mean(x) << endl;
    //cout << "y-mean = " << mean(y) << endl;

    //TF2 *circleFitFunction = new TF2("circleFitFunction","(x-[0])*(x-[0])+(y-[1])*(y-[1])-[2]*[2]",0,120,25,100);

    //TF1 *circleFitFunction = new TF1("circleFitFunction","sqrt([0]*[0]-(x-[1])*(x-[1]))+[2]",0,120);
    //circleFitFunction->GetXaxis()->SetRangeUser(0,xMax);
    //circleFitFunction->GetYaxis()->SetRangeUser(30,60);
    //circleFitFunction->SetParameters(-13,-134,2200);

    //cout << endl << "GRAPH JUST BEFORE FIT" << endl;
    //for(int i=0;i<g->GetN();i++)
    //    cout << g->GetX()[i] << " " << g->GetY()[i] << endl; //" " << g->GetZ()[i] << endl;

    cout << endl;

    //g->Fit("circleFitFunction","RW");

    //double x0=circleFitFunction->GetParameter(0);
    //double y0=circleFitFunction->GetParameter(1);
    //double r=circleFitFunction->GetParameter(2);
    //TEllipse *e = new TEllipse(x0,y0,r,r);
    //e->SetLineWidth(2);
    //e->SetFillStyle(0);
    //e->Draw();

    circle.Define(name2,x,y);

    //TCanvas *cCirc = new TCanvas();
    //cCirc->cd();
    //cout << "CURRENT CANVAS: " << gPad->GetCanvas()->GetName() << endl;
    //g->Draw("AL");
    //e->Draw();
    //cCirc->SaveAs("circleFit.C");
    //stringstream s;
    //s << "img/circles/step" << timestep << ".png";
    //cCirc->SaveAs(s.str().data());
    circle.Print();

    //Get chi2 scaled by number of points in fit
    return circle.GetChi2s();
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

//Split a string into a string vector of words
vector<string> strSplit(string str)
{
    int len=str.length();
    stringstream ss(str);
    int numWords=1;
    bool foundSpace=false;

    //Count number of words
    for(int ch=0;ch<len;ch++)
    {
        if(isspace(str[ch]))
            foundSpace=true;
        if(!isspace(str[ch])&&foundSpace)
        {
            numWords++;
            foundSpace=false;
        }
    }

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

    //Number of bins
    int nx=hist->GetNbinsX();

    //Scan columns from right to left
    for(int i=nx;i>0;i--)
    {
        if(hist->GetBinContent(i,j)>=0.5)
        {
            guess=hist->GetXaxis()->GetBinCenter(i);
            break;
        }
    }

    return guess;
}

//Given a TH1D and a two bin numbers, draw a line between the points and solve for where y=yc (y-0.5)
double solveLinear(TH1D *hist,int bin1,int bin2,double yc)
{
    //Find x & y values at bins
    double x1=hist->GetBinCenter(bin1);
    double x2=hist->GetBinCenter(bin2);
    double y1=hist->GetBinContent(bin1);
    double y2=hist->GetBinContent(bin2);

    //Equation of line
    double m = (y2-y1)/(x2-x1);
    double b = y1-m*x1;

    //Solve
    double xc = (yc-b)/m;
    return xc;
}

//Guess boundary location and width of water molecule for a single row or column
void guessTanhFit(TH1D* hist,double* fitBounds,double &ld,double &width,double &boundary,double &xlo,double &xhi)
{
    //f(x)=ld/2*(1-tanh(4*(x-x0)/w))
    //x=xlo=boundary+width/2 => yc=ld/2*(1-tanh(2))
    //x=xhi=boundary-width/2 => yc=ld/2*(1+tanh(2))
    //We will find these x values: xlo and xhi
    //From there, calculate width = hi-low
    //Since we are fitting ld, though, we don't yet know it.
    //We will use a value of 1, the ideal bulk density of water

    //Number of bins
    int n=hist->GetNbinsX();
    //Bin content (y value)
    double yc;

    //Assume ld=1
    ld=1;

    //Temporary variable for checking whether parameters are within limits
    double tmp;

    //Lower and upper edges for width
    double ylo=ld/2*(1-tanh(2));
    double yhi=ld/2*(1+tanh(2));

    //Whether edges and boundary have been found
    bool foundLower=false;
    bool foundUpper=false;
    bool foundBoundary=false;

    //Scan rows from top to bottom to find the boundary and width
    for(int i=n;i>0;i--)
    {
        yc=hist->GetBinContent(i);
        //Find the first bin with boundary density (0.5) then draw a line between this bin and the previous (above) and solve for where the line=0.5
        if( yc>0.5 && !foundBoundary )
        {
            tmp=solveLinear(hist,i,i+1,0.5);
            foundBoundary=true;
            if(fitBounds[4]<tmp && tmp<fitBounds[5])
                boundary=tmp;
        }

        //Look for lower edge
        if( yc>ylo && !foundLower )
        {
            xlo=solveLinear(hist,i,i+1,ylo);
            foundLower=true;
        }

        //Look for upper edge
        if( yc>yhi && !foundUpper )
        {
            xhi=solveLinear(hist,i,i+1,yhi);
            foundUpper=true;
        }
    }

    //Guess width if within bounds. Otherwise, revert to default guess (5)
    tmp=xlo-xhi;
    if(fitBounds[2]<tmp && tmp<fitBounds[3])
        width=tmp;
}

//Fit TH1D to tanh function, solve for x where f(x)=0.5
//Only take bins after and including startBin
//fitType should be "row", "col", or "mono"
double solveTanhFit(TH1D* hist,TF1* tanhFit,double* fitBounds,int startBin,int frameStep,double bulkEdge,char *fitType,double pos,TLegend* tanhLegend,TLine** tanhLines,TText **tanhTexts,TPaveText *tanhTextBox)
{
    double val;
    int n=hist->GetNbinsX();
    double *x = new double[n-startBin+1];
    double *y = new double[n-startBin+1];
    double startPoint=hist->GetBinCenter(startBin);
    bool draw=false;

    //Plot bounds
    double tanhRangeUser[4] = {0,120,0,3};

    //Initial guesses for width and boundary location in case they cannot be guessed for frame #1
    double ldGuess=1;
    double width=5;
    double boundary=bulkEdge;

    //Row or col?
    static int fitNum=0;
    static char *prevFitType=fitType;
    //Reset counter when changing types
    if(fitType!=prevFitType)
        fitNum=0;

    //cout << "solveTanhFit" << endl;
    //cout << endl;
    //cout << "n=" << n << endl;
    //cout << "startBin=" << startBin << endl;
    //cout << "startPoint=" << startPoint << endl;

    //Best guess based on counting, not fitting
    double lowGuess,hiGuess;
    guessTanhFit(hist,fitBounds,ldGuess,width,boundary,lowGuess,hiGuess);
    /*
    cout << endl;
    cout << "Guessed fit to be: " << endl;
    cout << "ldGuess=" << ldGuess << endl;
    cout << "width=" << width << endl;
    cout << "boundary=" << boundary << endl;
    cout << endl;
    */
    tanhFit->SetParameters(ldGuess,width,boundary);

    //Add hist values to vector, including only those after and including startBin
    for(int i=startBin;i<=n;i++)
    {
        x[i-startBin]=hist->GetBinCenter(i);
        y[i-startBin]=hist->GetBinContent(i);
        //cout << "Point " << i-startBin << ": (" << x[i-startBin] << "," << y[i-startBin] << ")" << endl;
    }
    //cout << "Points done!" << endl << endl;

    //Create graph
    TGraph* tanhPointsGraph = new TGraph(n-startBin,x,y);
    //cout << "tanhPointsGraph->GetN()=" << tanhPointsGraph->GetN() << endl;

    tanhPointsGraph->Fit(tanhFit,"Q");

    //Get Parameters
    double ld=tanhFit->GetParameter(0);
    double w=tanhFit->GetParameter(1);
    double c=tanhFit->GetParameter(2);

    /*
    cout << endl;
    cout << "Calculated fit to be: " << endl;
    cout << "ldGuess=" << ld << endl;
    cout << "width=" << w << endl;
    cout << "boundary=" << c << endl;
    cout << endl;
    */

    //Assume the fit is okay until proven otherwise
    //abs(c) => fit didn't blow up
    //ld>0.5 =>fit is valid and intersects 0.5=>column contains droplet
    //Otherwise => failure
    if(abs(c)<1000 && ld>0.5)
    {
        //Solve for x coordinate where the tanhFit(x)=0.5
        val=w/4*atanh(1-1.0/ld)+c;

        //If successful but too low, call it a failure
        if(val<=startPoint)
            val=-1;
    }
    else
    {
        //Histogram doesn't contain bulk water
        val=-1;
        //cout << "TanhFit failed for " << fitType << " " << fitNum << endl;
    }

    //if(frameStep==8810000)
    //draw=true;
    //else
    //    draw=false;

    //Save image
    if(draw && val>0)
    {
        //Canvas
        TCanvas *cTanh = new TCanvas();

        //Title
        stringstream title;

        //Draw hist
        hist->GetXaxis()->SetRangeUser(tanhRangeUser[0],tanhRangeUser[1]);
        hist->GetYaxis()->SetRangeUser(tanhRangeUser[2],tanhRangeUser[3]);
        hist->Draw();

        //Text - Where is this bin located?
        title.str("");
        title << fitType << " pos: " << pos;
        TText *posText = new TText(.45,.85,title.str().data());
        posText->SetNDC();
        posText->Draw();

        //Set annotation lines

        //Set lines to span plot box
        for(int i=1;i<7;i++)
        {
            //Vertical lines have identical x coordinates
            if(tanhLines[i]->GetX1()==tanhLines[i]->GetX2())
            {
                tanhLines[i]->SetY1(tanhRangeUser[2]);
                tanhLines[i]->SetY2(tanhRangeUser[3]);
            }
            //Horizontal lines do not
            else
            {
                tanhLines[i]->SetX1(tanhRangeUser[0]);
                tanhLines[i]->SetX2(tanhRangeUser[1]);
            }
        }

        //ld (horizontal)
        tanhLines[0]->SetY1(ld);
        tanhLines[0]->SetY2(ld);
        tanhLines[0]->Draw();
        //x0
        tanhLines[2]->SetX1(c);
        tanhLines[2]->SetX2(c);
        tanhLines[3]->Draw();
        //x0guess
        tanhLines[3]->SetX1(boundary);
        tanhLines[3]->SetX2(boundary);
        tanhLines[3]->Draw();
        //lowGuess
        tanhLines[4]->SetX1(lowGuess);
        tanhLines[4]->SetX2(lowGuess);
        tanhLines[4]->Draw();
        //hiGuess
        tanhLines[5]->SetX1(hiGuess);
        tanhLines[5]->SetX2(hiGuess);
        tanhLines[5]->Draw();
        //lowBin
        tanhLines[6]->SetX1(startPoint);
        tanhLines[6]->SetX2(startPoint);
        tanhLines[6]->Draw();
        //sol
        tanhLines[1]->SetX1(val);
        tanhLines[1]->SetX2(val);
        tanhLines[1]->Draw();

        //Set texts
        //pos
        title.str("");
        title << "pos: " << pos;
        tanhTexts[0]->SetText(0,0,title.str().data());
        //w_guess
        title.str("");
        title << "w_guess: " << width;
        tanhTexts[1]->SetText(0,0,title.str().data());
        //w
        title.str("");
        title << "w: " << w;
        tanhTexts[2]->SetText(0,0,title.str().data());
        //x0_guess
        title.str("");
        title << "x0_guess: " << boundary;
        tanhTexts[3]->SetText(0,0,title.str().data());
        //x0
        title.str("");
        title << "x0: " << c;
        tanhTexts[4]->SetText(0,0,title.str().data());
        title.str("");

        //Draw annotations
        tanhTextBox->Draw();
        tanhLegend->Draw();

        //Draw graph
        tanhPointsGraph->SetMarkerStyle(20);
        tanhPointsGraph->Draw("same p");

        /*
        //Draw solution
        TLine *solLine= new TLine(val,0,val,1);
        solLine->SetLineWidth(3);
        solLine->SetLineColor(kOrange);
        solLine->Draw();
        */

        //Title
        title.str("");
        title << "img/tanh/step" << frameStep << "_" << fitType << setw(3) << setfill('0') << fitNum << ".png";
        hist->SetTitle(title.str().data());

        //Save
        cTanh->SaveAs(title.str().data());

        delete cTanh;
    }

    //Update variables
    fitNum++;
    prevFitType=fitType;

    delete [] x;
    delete [] y;

    delete tanhPointsGraph;
    return val;
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

//Keep track of which atoms join the monolayer, and save their radius scaled by the base radius to the vector rScaledJoin
int monoFlux(vector<double> r,vector<double> z,double* monoLimits,double baseRadius,TH1D* rScaledJoin,TH1D* rScaledLeave,int &nMono)
{
    //Number of atoms in the simulation
    int numAtoms=z.size();

    //List of atoms currently in the monolayer
    vector<int> monoList;
    //List of atoms in the monolayer last timestep
    static    vector<int> prevMonoList;

    //Calculate flux since last step
    int flux = 0;

    //Count how many molecules should be in the monolayer based on the flux.
    static int expected = 0;

    //Save ids of atoms below the bulk-mono threshold
    for(int i=0;i<numAtoms;i++)
    {
        //if(i%(numAtoms/20)==0)
        //    cout << "#" << i << ": " << z[i] << endl;
        if((monoLimits[0]<=z[i]) && (z[i]<=monoLimits[1]))
            monoList.push_back(i);
    }

    //Number of atoms in the monolayer
    int numMonoAtoms=monoList.size();
    //Number of atoms in the monolayer last step
    int numPrevMonoAtoms=prevMonoList.size();

    //Save number of molecules in monolayer
    nMono = numMonoAtoms;


    cout << endl;
    cout << "Prev monolayer: " << numPrevMonoAtoms << endl;
    cout << "This monolayer: " << numMonoAtoms << endl;
    cout << endl;

    //For each atom in the monolayer, check whether it was there last step
    for(int i=0;i<numMonoAtoms;i++)
    {
        //Join
        //If it's new to the monolayer, save it's scaled radial position
        if( !isIn(monoList[i],prevMonoList) )
        {
            flux++;
            rScaledJoin->Fill(r[i]/baseRadius);
        }
    }

    //For each atom in the monolayer last time, check whether it's still there
    for(int i=0;i<numPrevMonoAtoms;i++)
    {
        //Leave
        //If has just left the monolayer, save it's scaled radial position
        if( !isIn(prevMonoList[i],monoList) )
        {
            flux--;
            rScaledLeave->Fill(r[i]/baseRadius);
        }
    }

    //Check whether the expected value matches the counted value
    /*
    expected += flux;
    cout << "Expected: " << expected;
    cout << " Counted: " << numMonoAtoms;
    if (expected == numMonoAtoms)
        cout << " OK!" << endl;
    else
        cout << " NOOOOOOOOOOOO!" << endl;
    */

    //Save monoList for next time
    prevMonoList=monoList;

    return flux;
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

//Round up to a multiple
int roundUp(int numToRound, int multiple)
{
        if (multiple == 0)
                    return numToRound;

            int remainder = numToRound % multiple;
                if (remainder == 0)
                            return numToRound;

                    return numToRound + multiple - remainder;
}

//Find z limits on monolayer
void findMonoLimits(TH1D *hWaterDens,double *monoLimits)
{
    int n = hWaterDens->GetNbinsX();
    bool foundPeak=false;
    bool foundDip=false;
    double dens=0;

    for(int i=1;i<=n;i++)
    {
        dens=hWaterDens->GetBinContent(i);

        //Find first peak above 1 - this is the beginning of the monolayer
        if(!foundPeak&&dens>1)
        {
            foundPeak=true;
            //cout << "Beginning monolayer at " << hWaterDens->GetBinLowEdge(i) << endl;
            monoLimits[0]=hWaterDens->GetBinLowEdge(i);
        }

        //if(foundPeak)
        //    cout << "dens=" << dens << " at " << hWaterDens->GetBinLowEdge(i) << endl;

        //if(foundPeak&&dens>1)
        //    cout << "Still in monolayer at " << hWaterDens->GetBinLowEdge(i) << endl;

        //Find first drop below 1 after peak - this is the end of the monolayer
        if(foundPeak&&dens<=1)
        {
            foundDip=true;
            //cout << "Monolayer ends at " << hWaterDens->GetBinLowEdge(i) << endl;
            monoLimits[1]=hWaterDens->GetBinLowEdge(i);
            break;
        }
    }

    if(foundDip)
        cout << "Found monolayer limits: " << monoLimits[0] << " " << monoLimits[1] << endl;
    else
        cout << "FAILED to locate monolayer." << endl;

}

