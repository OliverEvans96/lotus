//From "A Few Methods for Fitting Circles to Data"
//By Dale Umbach & Kerry N. Jones

//g++ `root-config --glibs --cflags` CircleFitClass.cpp circletest.cpp -o circletest.out && ./circletest.out

#include "CircleFitClass.h"

//Constructors
CircleFit::CircleFit(TH2D *givenHist)
{
    intersected=false;
    hist=givenHist;
    dataOut.open("CircleData.txt");
    dataOut << "x | y " << endl;
}

CircleFit::CircleFit(char* filename,TH2D* givenHist)
{
    ifstream inFile(filename);
    countLines(inFile);
    x.resize(n);
    y.resize(n);
    Fill(filename);
    inFile.close();
    intersected=false;
    hist=givenHist;
    x0=y0=0;
    r=1;
    Fit();

    /*
     cout << "double x[]={"<<x[0];
     for(int i=1;i<n;i++) cout << "," << x[i];
    cout << "};\ndouble y[]={"<<y[0];
     for(int i=1;i<n;i++) cout << "," << y[i];
        cout << "};" << endl;
        */
}

CircleFit::CircleFit(char* givenName,vector<double> xCoords,vector<double> yCoords)
{
    name=givenName;
    x=xCoords;
    y=yCoords;
    n=x.size();
    intersected=false;
    x0=y0=0;
    r=1;
//     for(int i=0;i<x.size();i++)
//         cout << "("<< x[i] << "," << y[i] << ")" << endl;
    Fit();
}

//Destructor
CircleFit::~CircleFit()
{
    //Save density histogram
    TCanvas *cCircDens = new TCanvas();
    cCircDens->cd();
    cCircDens->SaveAs("cCircDens.png");
    cCircDens->SaveAs("cCircDens.C");
    dataOut.flush();
    dataOut.close();
    delete cCircDens;
}

//Give name & points - same as constructor
void CircleFit::Define(char* givenName,vector<double> xCoords,vector<double> yCoords)
{
    name=givenName;
    x=xCoords;
    y=yCoords;
    n=x.size();
    cout << "CircleFit received " << n << " points" << endl;

    //Mirror points about y-axis if necessary
    //if(mirror)
    //{
        cout << "Mirroring points" << endl;
        x.resize(2*n);
        y.resize(2*n);
        for(int i=0;i<n;i++)
        {
            x[n+i]=-x[i];
            y[n+i]=y[i];
        }

        //Double number of points
        n*=2;
    //}

    intersected=false;
    x0=0;
    y0=0;
    r=1;
    cout << "About to fit" << endl;
    GuessFit();
}

//Get fitting error
double CircleFit::GetChi2s()
{
    return chi2s;
}

//Print the size of the arrays
void CircleFit::GetSize()
{
    cout << "n=" << n << endl;
    cout << "x.size()=" << x.size() << endl;
    cout << "y.size()=" << y.size() << endl;
}

//Fill coordinate vectors from file
void CircleFit::Fill(char* filename)
{
    ifstream inFile(filename);
    inFile.ignore(256,'\n');
    for(int i=0;i<n;i++)
    {
        inFile >> x[i] >> y[i];
        inFile.ignore(256,'\n');
    }

    inFile.close();
}

double CircleFit::LinearResidual(const double *X)
{
	// X = [m, b]
	sumsq = 0;
	// nLinFit is number of points to fit
	for(int i=0; i<xLinFit.size();i++)
		sumsq += pow(yLinFit[i] - X[0]*xLinFit[i] - X[1],2);

	return sumsq;
}

//Function to minimize
double CircleFit::SumOfSquares(const double *X)
{
    //X={x0,y0,r}
    sumsq=0;
    for(int i=0;i<n;i++)
    {
            //sumsq+=square(X[2]-sqrt(square(x[i]-X[0])+square(y[i]-X[1])));
            sumsq+=square(square(X[2])-square(x[i]-X[0])-square(y[i]-X[1]));
    }

    return sumsq;
}

//Fit circle to points with equal weights using MLS
void CircleFit::GuessFit()
{
    A=n*sum(mult(x,x))-square(sum(x));
    B=n*sum(mult(x,y))-sum(x)*sum(y);
    C=n*sum(mult(y,y))-square(sum(y));
    D=0.5*(n*sum(mult(x,mult(y,y)))-sum(x)*sum(mult(y,y))+n*sum(mult(x,mult(x,x)))-sum(x)*sum(mult(x,x)));
    E=0.5*(n*sum(mult(y,mult(x,x)))-sum(y)*sum(mult(x,x))+n*sum(mult(y,mult(y,y)))-sum(y)*sum(mult(y,y)));

    /*
    cout << "sum(x) = " << sum(x) << endl;
    cout << "sum(y) = " << sum(y) << endl;
    cout << "A = " << A << endl;
    cout << "B = " << B << endl;
    cout << "C = " << C << endl;
    cout << "D = " << D << endl;
    cout << "E = " << E << endl;
    */

    //Center
    x0=(D*C-B*E)/(A*C-B*B);
    y0=(A*E-B*D)/(A*C-B*B);

    //Radius
    findRadius();

    //Report
    cout << "Using MLS: (" << setprecision(10) << x0 << "," << y0 << "," << r << ")" << endl;
}

//Fit circle to points using numerical minimization
void CircleFit::Fit()
{

    cout << "Points:" << endl;
     for(int i=0;i<x.size();i++)
         cout << "("<< x[i] << "," << y[i] << ")" << endl;

    //Initialize minimizer
    minimizer.SetMaxFunctionCalls((int)1e7);
    minimizer.SetMaxIterations((int)1e7);
    minimizer.SetTolerance(1e-8);

    //Step size & initial values for minimization
    stepVal=1e-3;
    step[0]=stepVal;
    step[1]=stepVal;
    step[2]=stepVal;
    init[0]=x0;
    init[1]=y0;
    init[2]=r;

    //Wrap function
    ROOT::Math::Functor functor(this,&CircleFit::SumOfSquares,3);
    minimizer.SetFunction(functor);

    //Initialize variables
    //cout << "init" << endl;
    minimizer.SetLimitedVariable(0,"x0",init[0],step[0],-200,20);
    minimizer.SetLimitedVariable(1,"y0",init[1],step[1],-200,20);
    minimizer.SetLimitedVariable(2,"r",init[2],step[2],10,250);

    //Calculate center without weights
    GuessFit();

    //Minimize!
    minimizer.Minimize();

    //Get values
    x0=minimizer.X()[0];
    y0=minimizer.X()[1];
    r=abs(minimizer.X()[2]);

    //Save chi2 (scaled by number of points
    chi2s=minimizer.MinValue()/(r*r*n);
    //cout << "chi2s:" << chi2s << endl;

    //Report
    //cout << "Using minimizer: (" << setprecision(10) << x0 << "," << y0 << "," << r << ")" << endl;
}

//Find the largest x value at which the circle intersects a constant c
double CircleFit::Intersect(double c)
{
    y1=c;
    intersected=true;
    x1=sqrt(square(r)-square(c-y0))+x0;
    return x1;

}

//Calculate contact angle
double CircleFit::ContactAngle()
{
    if(intersected)
    {
        cosTheta=(y1-y0)/r;
        return cosTheta;
    }

    else
        cout << "Call Intersect(double c) before ContactAngle()" << endl;
}

//Contact angle from points within 5A of the rightmost point.
double CircleFit::LinearContactAngle()
{
	// Cutoff from rightmost point (5A)
	double cutoff = 5;
	// Points within cutoff
	vector<double> xLinFit, yLinFit;
	// Assume that points are ordered by x value
	double xmax = x[n-1];
	if(intersected)
	{
		// Loop through points, checking cutoff
		for(int i=0; i<n; i++)
		{
			if(abs(x[i] - xmax) < cutoff)
			{
				xLinFit.push_back(x[i]);
				yLinFit.push_back(y[i]);
			}
		}
	}

    //Initialize minimizer
    linMin.SetMaxFunctionCalls((int)1e7);
    linMin.SetMaxIterations((int)1e7);
    linMin.SetTolerance(1e-8);

    //Wrap function
    ROOT::Math::Functor functor(this,&CircleFit::LinearResidual,3);
    linMin.SetFunction(functor);

    //Initialize variables
    //cout << "init" << endl;
	double m_init = -1;
	double b_init = CircleFit::Height();
    linMin.SetLimitedVariable(0,"m",init[0],step[0],-10,10);
    linMin.SetLimitedVariable(1,"b",init[1],step[1],-200,200);

    //Minimize!
    linMin.Minimize();

    //Get values
    m_lin = linMin.X()[0];
    b_lin = linMin.X()[1];
}

//Draw tangent line
TLine* CircleFit::DrawTangentLine()
{
    if(intersected)
    {
        //Calculate slope & y-intercept
        //Negative reciprocal of radius slope
        m=(x0-x1)/(y1-y0);
        b=y1-m*x1;

        //Determine    which boundaries to use
        //Line will be drawn from (x2,y2) to (x3,y3)

        //x2 on bottom
        if(m*xlo+b<ylo)
            x2=(ylo-b)/m;
        //x2 on top
        else if(m*xlo+b>yhi)
            x2=(yhi-b)/m;
        //x2 on left
        else
            x2=xlo;

        //y2 is easy
        y2=m*x2+b;

        //x3 on bottom
        if(m*xhi+b<ylo)
            x3=(ylo-b)/m;
        //x3 on top
        else if(m*xhi+b>yhi)
            x3=(yhi-b)/m;
        //x3 on right
        else
            x3=xhi;

        //y3 is easy
        y3=m*x3+b;

        TLine *tangentLine = new TLine(x2,y2,x3,y3);
        tangentLine->SetLineColor(kViolet);
        tangentLine->SetLineWidth(3);
        tangentLine->Draw();

        return tangentLine;
    }
    else
        cout << "Call Intersect(double c) before ContactAngle()" << endl;
}

//Calculate the height of the circle at x=0
double CircleFit::Height()
{
    //Take the best of two definitions of height
    double h1,h2;

    //Intersection w/ y-axis
    h1=sqrt(square(r)-square(x0))+y0;
    //Highest point on circle
    //Don't use this one if the center is left of the y-axis
    if(x0>0)
        h2=y0+r;
    else
        h2=0;

    //Return maximum
    if(h1<h2)
        height=h2;
    else
        height=h1;

    return height;
}

//Draw points and fitted circle
TEllipse* CircleFit::Draw(bool drawPoints)
{
    //PI
    PI = 3.141592653589793;

    //Get margins
    margins[0]=gPad->GetLeftMargin();
    margins[1]=gPad->GetBottomMargin();
    margins[2]=gPad->GetRightMargin();
    margins[3]=gPad->GetTopMargin();

    //Get limits (these include values in the margins);
    //gPad->GetCanvas()->GetRange(xlo,ylo,xhi,yhi);

    /*
    cout << endl;
    cout << "Current canvas: " << gPad->GetCanvas()->GetName() << endl;
    cout << "Canvas range 1: " << xlo << " " << ylo << " " << xhi << " " << yhi << endl;

    //Account for margins
    xlo += (xhi-xlo)*margins[0];
    ylo += (xhi-xlo)*margins[1];
    xhi -= (yhi-ylo)*margins[2];
    yhi -= (yhi-ylo)*margins[3];
    */
    xlo=hist->GetXaxis()->GetXmin();
    ylo=hist->GetYaxis()->GetXmin();
    xhi=hist->GetXaxis()->GetXmax();
    yhi=hist->GetYaxis()->GetXmax();

    //Copy vectors to arrays
    double* xA=&x[0];
    double* yA=&y[0];

    //Points
    if(drawPoints)
    {
        TGraph* g = new TGraph(n,xA,yA);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1);
        g->Draw("AP");
    }

    //minAngle=0;
    //maxAngle=360;
    /*
    //Find min and max angle to draw
    //If circle intersects bottom
    if(ylo-y0<r)
        minAngle=asin((ylo-y0)/r)*180/PI;
    //If circle intersects left
    if(xlo-x0<r)
        maxAngle=acos((xlo-x0)/r)*180/PI;
    */

    //Fitted Circle
    TEllipse* e = new TEllipse(x0,y0,r,r);
    e->SetLineWidth(3);
    e->SetFillStyle(0);
    //e->SetNoEdges();
    e->Draw("same");

    return e;
//     stringstream fname;
//     fname << "img/circles/" << name << ".jpg";

//     c1->SaveAs(fname.str().data());
//     delete c1,g,e;
}

double CircleFit::sum(vector<double> v)
{
    double sum=0;
    int n=v.size();
    for(int i=0;i<n;i++)
        sum+=v[i];
    return sum;
}

//Square
double CircleFit::square(double x)
{
    return x*x;
}

//Elementwise vector addition
vector<double> CircleFit::add(vector<double> u,vector<double> v)
{
    vector<double> w=u;
    int n=v.size();
    for(int i=0;i<n;i++)
        w[i]+=v[i];
    return w;
}

//Elementwise vector multiplication
vector<double> CircleFit::mult(vector<double> u,vector<double> v)
{
    vector<double> w=u;
    int n=v.size();
    for(int i=0;i<n;i++)
        w[i]*=v[i];
    return w;
}

//Calculate mean
double CircleFit::mean(vector<double> x)
{
    double sum=0;
    for(int i=0;i<x.size();i++)
        sum+=x[i];
    return sum/x.size();
}

//Calculate standard deviation
double CircleFit::stddev(vector<double> x)
{
    double m = mean(x);
    double sum=0;
    for(int i=0;i<x.size();i++)
        sum+=(m-x[i])*(m-x[i]);
    return sqrt(sum/(x.size()-1));
}

//Arctanh
double CircleFit::atanh(double x)
{
    return log((1+x)/(1-x))/2;
}

//Find radius from points and center
void CircleFit::findRadius()
{
    int n=x.size();
    r=0;

    for(int i=0;i<n;i++)
        r+=sqrt(square(x[i]-x0)+square(y[i]-y0))/n;
}

///Get x center
double CircleFit::GetXCenter()
{
    return x0;
}

//Get y center
double CircleFit::GetYCenter()
{
    return y0;
}

//Get radius
double CircleFit::GetRadius()
{
    return r;
}

//Print x & y centers & radius
void CircleFit::Print()
{
    cout << "X Center: " << x0 << endl;
    cout << "Y Center: " << y0 << endl;
    cout << "Radius: " << r << endl;
    cout << endl;
}

//Count number of lines in file
void CircleFit::countLines(ifstream &inFile)
{
    //Beginning of file
    inFile.seekg(0,ios::beg);

    //Count number of atoms
    bool countFlag=true;
    string line;
    int numLines=0;

    //Ignore the first line
    getline(inFile,line);

    while(countFlag)
    {
        getline(inFile,line);

        if(inFile.eof())
            countFlag=false;
        else
            numLines++;
    }

    //Unset eof flag (if set) & return to beginning of file
    inFile.clear();
    inFile.seekg(0,ios::beg);

    n=numLines-1;
}

//Check whether point is already in graph
bool CircleFit::inGraph(TGraph *g,double xCheck,double yCheck)
{
    bool in=false;
    double *xGraph=g->GetX();
    double *yGraph=g->GetY();

    for(int i=0;i<g->GetN();i++)
    {
        if( (xCheck==xGraph[i]) && (yCheck==yGraph[i]) )
            in=true;
    }

    return in;
}
