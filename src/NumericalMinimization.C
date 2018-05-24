#include "Math/GSLSimAnMinimizer.h"
#include "Math/Functor.h"
 
double RosenBrock(const double *xx )
{
	const Double_t x = xx[0];
	const Double_t y = xx[1];
	return -exp(-(x-3)*(x-3))-exp(-(y+2)*(y+2));
}

 
int NumericalMinimization()
{
	ROOT::Math::GSLSimAnMinimizer min;
 
	min.SetMaxFunctionCalls(1000000);
	min.SetMaxIterations(100000);
	min.SetTolerance(1e-7);
 
	double step[2] = {0.01,0.01};
	double variable[2] = { -1,1.2};
 
	ROOT::Math::Functor f(&RosenBrock,2); 

	min.SetFunction(f);
 
	// Set the free variables to be minimized!
	min.SetVariable(0,"x",variable[0], step[0]);
	min.SetVariable(1,"y",variable[1], step[1]);
 
	min.Minimize(); 
 
	const double *xs = min.X();
	cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " 
			<< RosenBrock(xs) << endl;
 
	return 0;
}
