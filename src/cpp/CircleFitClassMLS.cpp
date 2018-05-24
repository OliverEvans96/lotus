//From "A Few Methods for Fitting Circles to Data"
//By Dale Umbach & Kerry N. Jones

//g++ `root-config --glibs --cflags` CircleFitClass.cpp circletest.cpp -o circletest.out && ./circletest.out

#include "CircleFitClassMLS.h"


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
	x0=y0=0;
	r=1;
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
	x0=y0=0;
	r=1;
// 	for(int i=0;i<x.size();i++)
// 		cout << "("<< x[i] << "," << y[i] << ")" << endl;
	Fit();
}

//Destructor
CircleFit::~CircleFit()
{
}

//Give name & points - same as constructor
void CircleFit::Define(char* givenName,vector<double> xCoords,vector<double> yCoords)
{
	name=givenName;
	x=xCoords;
	y=yCoords;
	n=x.size();
	intersected=false;
	x0=y0=0;
	r=1;
// 	for(int i=0;i<x.size();i++)
// 		cout << "("<< x[i] << "," << y[i] << ")" << endl;
	Fit();
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
	A=n*sum(mult(x,x))-square(sum(x));
	B=n*sum(mult(x,y))-sum(x)*sum(y);
	C=n*sum(mult(y,y))-square(sum(y));
	D=0.5*(n*sum(mult(x,mult(y,y)))-sum(x)*sum(mult(y,y))+n*sum(mult(x,mult(x,x)))-sum(x)*sum(mult(x,x)));
	E=0.5*(n*sum(mult(y,mult(x,x)))-sum(y)*sum(mult(x,x))+n*sum(mult(y,mult(y,y)))-sum(y)*sum(mult(y,y)));
	
	//Center
	tmpx=(D*C-B*E)/(A*C-B*B);
	tmpy=(A*E-B*D)/(A*C-B*B);

	x0=tmpx;
	y0=tmpy;

	//Radius
	findRadius();
	r=tmpr;

	/*
	//Verify
	if( (abs(tmpx)<=25) && (abs(tmpy)<=25) && (tmpr<=100) )
	{
		x0=tmpx;
		y0=tmpy;
		r=tmpr;
		cout << endl;
		cout << "Circle Fit SUCCESS" << endl;
		cout << "x0=" << x0 << endl;
		cout << "y0=" << y0 << endl;
		cout << "r=" << r << endl;
		cout << endl;
	}
	else
	{
		cout << "Circle Fit FAILURE" << endl;
		cout << "x0=" << x0 << endl;
		cout << "y0=" << y0 << endl;
		cout << "r=" << r << endl;
		cout << endl;
	}
	*/
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
		cosTheta=(x1-x0)/r;
		return cosTheta;
	}
	
	else
		cout << "Call Intersect(double c) before ContactAngle()" << endl;
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

		//Determine	which boundaries to use
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
		tangentLine->SetLineWidth(2);
		tangentLine->Draw();

		return tangentLine;
	}
	else
		cout << "Call Intersect(double c) before ContactAngle()" << endl;
}

//Calculate the height of the circle at x=0
double CircleFit::Height()
{
	height=sqrt(square(r)-square(x0))+y0;
	return height;
}

//Draw points and fitted circle
TEllipse* CircleFit::Draw(bool drawPoints)
{
	//PI
	PI = 3.141592653589793;

	//Get margins
	margins[0]=(double) gPad->GetLeftMargin();
	margins[1]=(double) gPad->GetBottomMargin();
	margins[2]=(double) gPad->GetRightMargin();
	margins[3]=gPad->GetTopMargin();

	//Get limits (these include values in the margins);
	gPad->GetRange(xlo,ylo,xhi,yhi);

	//Account for margins
	xlo += (xhi-xlo)*margins[0];
	ylo += (xhi-xlo)*margins[1];
	xhi -= (yhi-ylo)*margins[2];
	yhi -= (yhi-ylo)*margins[3];

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
	
	//Find min and max angle to draw
	//If circle intersects bottom
	if(ylo-y0<r)
		minAngle=asin((ylo-y0)/r)*180/PI;
	else
		minAngle=0;
	//If circle intersects left
	if(xlo-x0<r)
		maxAngle=acos((xlo-x0)/r)*180/PI;
	else
		maxAngle=360;

	//Fitted Circle
	TEllipse* e = new TEllipse(x0,y0,r,r,minAngle,maxAngle);
	e->SetLineWidth(3);
	e->SetFillStyle(0);
	e->SetNoEdges();
	e->Draw("same");

	return e;
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
	tmpr=0;
	
	for(int i=0;i<n;i++)
		tmpr+=sqrt(square(x[i]-x0)+square(y[i]-y0))/n;
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

