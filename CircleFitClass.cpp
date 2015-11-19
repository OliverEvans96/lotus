//From "A Few Methods for Fitting Circles to Data"
//By Dale Umbach & Kerry N. Jones

//g++ `root-config --glibs --cflags` CircleFitClass.cpp circletest.cpp -o circletest.out && ./circletest.out

#include "CircleFitClass.h"

//Constructors
CircleFit::CircleFit()
{
	intersected=false;
}

CircleFit::CircleFit(char* filename)
{
	ifstream inFile(filename);
	countLines(inFile);
	x.resize(n);
	y.resize(n);
	Fill(filename);
	inFile.close();
	intersected=false;
	Fit();
	
// 	cout << "double x[]={"<<x[0];
// 	for(int i=1;i<n;i++) cout << "," << x[i];
// 	cout << "};\ndouble y[]={"<<y[0];
// 	for(int i=1;i<n;i++) cout << "," << y[i];
// 	cout << "};" << endl;
}

CircleFit::CircleFit(char* givenName,vector<double> xCoords,vector<double> yCoords)
{
	name=givenName;
	x=xCoords;
	y=yCoords;
	n=x.size();
	intersected=false;
// 	for(int i=0;i<x.size();i++)
// 		cout << "("<< x[i] << "," << y[i] << ")" << endl;
	Fit();
}

//Destructor
CircleFit::~CircleFit()
{
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

//Fit circle to points
void CircleFit::Fit()
{
	double A=n*sum(mult(x,x))-square(sum(x));
	double B=n*sum(mult(x,y))-sum(x)*sum(y);
	double C=n*sum(mult(y,y))-square(sum(y));
	double D=0.5*(n*sum(mult(x,mult(y,y)))-sum(x)*sum(mult(y,y))+n*sum(mult(x,mult(x,x)))-sum(x)*sum(mult(x,x)));
	double E=0.5*(n*sum(mult(y,mult(x,x)))-sum(y)*sum(mult(x,x))+n*sum(mult(y,mult(y,y)))-sum(y)*sum(mult(y,y)));
	
	//Center
	x0=(D*C-B*E)/(A*C-B*B);
	y0=(A*E-B*D)/(A*C-B*B);
	
	//Radius
	findRadius();
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
		theta=atanh( (x1-x0)/(y1-y0) );
		cosTheta=(x1-x0)/r;
		return cosTheta;
	}
	
	else cout << "Call Intersect(double c) before ContactAngle()" << endl;
}

//Calculate the height of the circle at x=0
double CircleFit::Height()
{
	height=sqrt(square(r)-square(x0))+y0;
	return height;
}

//Draw points and fitted circle
void CircleFit::Draw(bool drawPoints)
{
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
	
	//Fitted Circle
	TEllipse* e = new TEllipse(x0,y0,r,r);
	e->SetLineWidth(3);
	e->SetFillStyle(0);
	e->Draw("same");
	
// 	stringstream fname;
// 	fname << "img/circles/" << name << ".jpg";
	
// 	c1->SaveAs(fname.str().data());
// 	delete c1,g,e;
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

