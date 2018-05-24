#include "RConfigure.h"
#ifdef R__HAS_MATHMORE
#include "Math/MultiRootFinder.h"
#endif
#include "Math/WrappedMultiTF1.h"
#include "TF2.h"
#include "TError.h"
// example of using multi root finder based on GSL
// need to use an algorithm not requiring the derivative
//like hybrids (default), hybrid, dnewton, broyden
using namespace ROOT::Math;
void exampleMultiRoot(const char * algo = 0, int printlevel = 1) {
#ifndef R__HAS_MATHMORE
   Error("exampleMultiRoot","libMathMore is not available - cannot run this tutorial");
#else
   ROOT::Math::MultiRootFinder r(algo);
    //defining the function
   // use Rosenbrock functions
   TCanvas *c = new TCanvas();
   TF2 * f1 = new TF2("f1","[0]*(1-x)+[1]*y");
   f1->SetRange(-10,10);
   f1->Draw();
   TF2 * f2 = new TF2("f2","[0]*(y-x*x)");
   f2->SetRange(-10,10);
   f2->Draw();

   f1->SetParameters(1,0);
   f2->SetParameter(0,10);
   // wrap the functions
   ROOT::Math::WrappedMultiTF1 g1(*f1,2);
   ROOT::Math::WrappedMultiTF1 g2(*f2,2);
   r.AddFunction(g1);
   r.AddFunction(g2);
   r.SetPrintLevel(printlevel);
   // starting point
   double x0[2]={-1,-1};
   r.Solve(x0);
#endif
}

