#include "CircleFitClass.h"
#include <fstream>
#include <iostream>

//g++ `root-config --glibs --cflags` CircleFitClass.cpp circletest.cpp -o circletest.out && ./circletest.out

using namespace std;

int main()
{
	double x[]={1.0,7.0,6.0,7.0};
	double y[]={9.0,4.0,7.0,2.0};
	
// 	double x[]={94.4602,92.4716,84.517,80.5398,72.5852,64.6307,58.6648,54.6875,46.733,36.7898};
// 	double y[]={30.9605,32.9342,34.9079,36.8816,38.8553,40.8289,42.8026,44.7763,46.75,48.7237};
	
	vector<double> xV(x,x+4);
	vector<double> yV(y,y+4);
	
	CircleFit Test(xV,yV);
// 	CircleFit Test("points.txt");
	
	Test.GetSize();
	Test.Print();
	
// 	cout << Test.Intersect(3) << endl;
	Test.Draw();

	return 0;
}
