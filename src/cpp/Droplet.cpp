#include "Droplet.h"

////////////////////////
// Droplet Components //
////////////////////////

Monolayer::Monolayer() {}
Monolayer::~Monolayer() {
  // TODO: Pretty weird to allocate in Droplet
  // and delete in Monolayer.
  delete hMono;
}

void Monolayer::setContext(Options _options, SimData *_simDataPtr, AtomArray *_atomArrayPtr) {
  options = _options;
  simDataPtr = _simDataPtr;
  atomArrayPtr = _atomArrayPtr;
}

void Monolayer::fillOne(Atom &atom) {
  int mass;
  mass = simDataPtr->masses[atom.type];
  // TODO: Make this more efficient by calculating
  // only the necessary "radius" (y or r).
  if(options.geometry == "spherical") {
    hMono->Fill(atom.r, mass);
    //cout << "fill mono @ " << atom.r << endl;
  }
  else if(options.geometry == "cylindrical") {
    hMono->Fill(atom.y, mass);
    //cout << "fill mono @ " << atom.r << endl;
  }
}

void Monolayer::fill(AtomArray &atoms) {
  Atom atom;
  // Add to monolayer
  for(int i=0; i<simDataPtr->numAtoms; i++) {
    atoms.getAtom(i, atom);
    if(isIn(atoms.type[i], simDataPtr->liquidTypes)) {
      if(inMonolayer(atom)) {
        fillOne(atom);
      }
    }
  }
}

void Monolayer::convertUnits() {
  // Must be called after monoLimits are found.
  double dv;

  // Divide by number of steps per frame
  hMono->Scale(1.0/simDataPtr->stepsPerFrame);
  // Calculate monolayer bin volume
  // options.dv is for hDroplet, which uses options.dz.
  dv = options.dv * (zlim[1] - zlim[0])/options.dz;
  // Divide by volume to get density
  hMono->Scale(1.0/dv);
  // Convert units from amu/AA to g/cc
  hMono->Scale(NANO_DENS_TO_MACRO);
}

// Whether an atom is in the monolayer
bool Monolayer::inMonolayer(Atom &atom) {
  return ((zlim[0] < atom.z) && (atom.z > zlim[1]));
}

//Keep track of which atoms join the monolayer, and save their radius scaled by the base radius to the vector rScaledJoin
int Monolayer::monoFlux(vector<double> r,vector<double> z,double* monoLimits,double baseRadius,TH1D* rScaledJoin,TH1D* rScaledLeave,int &nMono)
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

//Find z limits on monolayer
void Monolayer::findMonoLimits(TH1D *hLiquidDens,double *monoLimits)
{
    int n = hLiquidDens->GetNbinsX();
    bool foundPeak=false;
    bool foundDip=false;
    double dens=0;

    for(int i=1;i<=n;i++)
    {
        dens=hLiquidDens->GetBinContent(i);

        //Find first peak above 1 - this is the beginning of the monolayer
        if(!foundPeak&&dens>1)
        {
            foundPeak=true;
            //cout << "Beginning monolayer at " << hLiquidDens->GetBinLowEdge(i) << endl;
            monoLimits[0]=hLiquidDens->GetBinLowEdge(i);
        }

        //if(foundPeak)
        //    cout << "dens=" << dens << " at " << hLiquidDens->GetBinLowEdge(i) << endl;

        //if(foundPeak&&dens>1)
        //    cout << "Still in monolayer at " << hLiquidDens->GetBinLowEdge(i) << endl;

        //Find first drop below 1 after peak - this is the end of the monolayer
        if(foundPeak&&dens<=1)
        {
            foundDip=true;
            //cout << "Monolayer ends at " << hLiquidDens->GetBinLowEdge(i) << endl;
            monoLimits[1]=hLiquidDens->GetBinLowEdge(i);
            break;
        }
    }

    if(foundDip)
        cout << "Found monolayer limits: " << monoLimits[0] << " " << monoLimits[1] << endl;
    else
        cout << "FAILED to locate monolayer." << endl;

}

CircularBulk::CircularBulk() {
  gCirclePoints = new TGraph();
}

CircularBulk::~CircularBulk() {
  delete gCirclePoints;
}

void CircularBulk::setContext(Options _options, SimData *_simDataPtr, AtomArray *_atomArrayPtr) {
  options = _options;
  simDataPtr = _simDataPtr;
  atomArrayPtr = _atomArrayPtr;
}

//Guess boundary of water molecule by counting for a single row
double CircularBulk::guessRowBoundary(TH2D* hist,int j)
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
//Find the edge of droplet by tanh fitting for each row given TH2D
void CircularBulk::findBoundaryPoints(TH2D* hist,TGraph *circlePointsGraph,char* aOrR,double *monoLimits,TH1D *hMono,double& rBulkMax,double &monoEdge,double &bulkEdge,int frameStep,TLegend* tanhLegend,TLine** tanhLines,TText **tanhTexts,TPaveText *tanhTextBox)
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
    // TODO: Set function limits from options, not hard coded
    TF1* tanhFit = new TF1("tanhFit","[0]/2*(1-tanh(4*(x-[2])/([1])))",0,300);

    //Set Bounds on parameters
    // TODO: Set bounds from options
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
        cout << "RF j = " << j << endl;
        //Create projection
        px = hist->ProjectionX("px",j,j);
        center = hist->GetYaxis()->GetBinCenter(j);

        //Solve for x coordinate where tanhFit(x)=0.5
        //Include all bins
        tmpx=solveTanhFit(px,tanhFit,fitBounds,1,frameStep,bulkEdge,"row",center,tanhLegend,tanhLines,tanhTexts,tanhTextBox);
        cout << "tmpx = " << tmpx << endl;
        //cout << "Row " << j << ": tmpx=" << tmpx << endl;
        if(tmpx>0)
        {
            //Found another row containing droplet
            tmpy=hist->GetYaxis()->GetBinCenter(j);
            circlePointsGraph->SetPoint(n,tmpx,tmpy);
            printf("CPG Set %d @ (%.2f, %.2f).\n", n, tmpx, tmpy);

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

    delete tanhFit;
}

//Given a graph of points, a maximum x value, and the y coordinate of the interface with the substrate, fit a circle to the points and find the intersection of the circle with the substrate interface. The result is the bulk-monolayer interface
double CircularBulk::fitCircle(TGraph* gCirclePoints,CircleFit &circle,double xMax,int timestep) {
    //gCirclePoints->Draw("SAME");
    //Fit circle and intersect it with a constant c (Interface)
    char* name1 = (char*) gCirclePoints->GetTitle();
    stringstream nameStream;
    nameStream << name1 << setw(8) << setfill('0') << timestep;
    char* name2 = (char*) nameStream.str().data();
    int n=gCirclePoints->GetN();
    vector<double> x(n),y(n);
    double xTest,yTest;

    cout << "xMax=" << xMax << endl;

    //Sort points
    gCirclePoints->Sort();

//     //Limits for circle fitting for first 50 timesteps
//     double lowLim=30;
//     double highLim=50;
//
    //Number of points with x<xMax
    int nValid=0;

    for(int i=0;i<n;i++)
    {
        //Test the point before including it
        //gCirclePoints->GetPoint(i,xTest,yTest);
        xTest=gCirclePoints->GetX()[i];
        yTest=gCirclePoints->GetY()[i];

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

    cout << endl;

    circle.Define(name2,x,y);

    circle.Print();

    //Get chi2 scaled by number of points in fit
    return circle.GetChi2s();
}

//Guess boundary location and width of water molecule for a single row or column
void CircularBulk::guessTanhFit(TH1D* hist,double* fitBounds,double &ld,double &width,double &boundary,double &xlo,double &xhi)
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
double CircularBulk::solveTanhFit(TH1D* hist, TF1* tanhFit, double* fitBounds, int startBin, int frameStep, double bulkEdge, string fitType, double pos, TLegend* tanhLegend, TLine** tanhLines, TText **tanhTexts, TPaveText *tanhTextBox)
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
    static string prevFitType=fitType;
    //Reset counter when changing types
    if(fitType!=prevFitType)
        fitNum=0;

    //cout << "solveTanhFit" << endl;
    //cout << endl;
    //cout << "n=" << n << endl;
    //cout << "startBin=" << startBin << endl;
    //cout << "startPoint=" << startPoint << endl;

    cout << "THF" << endl;
    //Best guess based on counting, not fitting
    double lowGuess,hiGuess;
    guessTanhFit(hist,fitBounds,ldGuess,width,boundary,lowGuess,hiGuess);
    cout << "Guess:" << endl;
    cout << "ldGuess = " << ldGuess << endl;
    cout << "width = " << width << endl;
    cout << "boundary = " << boundary << endl;
    cout << "lowGuess = " << lowGuess << endl;
    cout << "hiGuess = " << hiGuess << endl;
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
        if(y[i-startBin] != 0) {
          cout << "(nonzero Point) " << i-startBin << ": (" << x[i-startBin] << "," << y[i-startBin] << ")" << endl;
        }
    }
    //cout << "Points done!" << endl << endl;

    //Create graph
    TGraph* tanhPointsGraph = new TGraph(n-startBin,x,y);
    //cout << "tanhPointsGraph->GetN()=" << tanhPointsGraph->GetN() << endl;

    tanhPointsGraph->Fit(tanhFit,"Q");

    //Get Parameters
    double ld=tanhFit->GetParameter(0);
    double w=tanhFit->GetParameter(1);
    double x0=tanhFit->GetParameter(2);

    cout << "Fit Params:" << endl;
    cout << "ld = " << ld << endl;
    cout << "w = " << w << endl;
    cout << "x0 = " << x0 << endl;

    /*
    cout << endl;
    cout << "Calculated fit to be: " << endl;
    cout << "ldGuess=" << ld << endl;
    cout << "width=" << w << endl;
    cout << "boundary=" << x0 << endl;
    cout << endl;
    */

    //ld>0.5 =>fit is valid and intersects 0.5=>column contains droplet
    //Otherwise => failure
    // TODO: Set max in options
    if(abs(x0)<1000 && ld>0.0)
    {
      // Use point where dens = ld/2.
      val = x0;

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

    cout << "val = " << val << endl;

    //Update variables
    fitNum++;
    prevFitType=fitType;

    delete [] x;
    delete [] y;

    delete tanhPointsGraph;
    return val;
}

Droplet::Droplet(AtomArray &atomArray) {
  // TODO: Determine radius from options?
  rDensCyl = 20;

  setContext(atomArray);
  createHists();
}

Droplet::~Droplet() {
  delete hDroplet;
  delete cDroplet;
}

void Droplet::setContext(AtomArray &atomArray) {
  atomArrayPtr = &atomArray;
  simDataPtr = atomArrayPtr->simDataPtr;
  options = simDataPtr->options;

  bulk.setContext(options, simDataPtr, atomArrayPtr);
  monolayer.setContext(options, simDataPtr, atomArrayPtr);
}

void Droplet::fillOne(Atom &atom) {
  double mass, z;

  mass = simDataPtr->masses[atom.type];
  // z relative to top of substrate
  z = atom.z - simDataPtr->substrateTop;
  // cout << "atom.z = " << atom.z << endl;
  // cout << "substrateTop = " << simDataPtr->substrateTop << endl;
  // cout << "z = " << z << endl;

  // 2D Density
  hDroplet->Fill(atom.r, z, mass);
  // 1D Density
  if(atom.r < rDensCyl) {
    hLiquidDens->Fill(atom.z, mass);
  }
}

void Droplet::fill(AtomArray &atoms) {
  Atom atom;
  reset();
  cout << "Filling w.r.t. z=" << simDataPtr->substrateTop << endl;
  cout << "monoLimits: " << monolayer.zlim[0] << " " << monolayer.zlim[1] << endl;
  for(int i=0; i<simDataPtr->numAtoms; i++) {
    // If liquid
    if(isIn(atoms.type[i], simDataPtr->liquidTypes)) {
      atoms.getAtom(i, atom);
      atom.calculateNonCartesian();
      fillOne(atom);
    }
  }

  convertUnits();
  monolayer.convertUnits();
}

void Droplet::findMonolayer() {
  monolayer.findMonoLimits(hLiquidDens, monolayer.zlim);

  cout << "ML: hLD @" << hLiquidDens << endl;
}

void Droplet::reset() {
  hDroplet->Reset();
  hLiquidDens->Reset();
}

void Droplet::convertUnits() {
  // Divide by number of steps per frame
  hDroplet->Scale(1.0/simDataPtr->stepsPerFrame);
  hLiquidDens->Scale(1.0/simDataPtr->stepsPerFrame);
  // Divide by volume to get density
  hDroplet->Scale(1.0/options.dv);
  hLiquidDens->Scale(1.0/(PI*options.dz*square(rDensCyl)));
  // Convert units from amu/AA^3 to g/cc
  hDroplet->Scale(NANO_DENS_TO_MACRO);
  hLiquidDens->Scale(NANO_DENS_TO_MACRO);
}

void Droplet::createHists() {
  double zlo, zhi;
  double rlo, rhi;
  // Set limits on histograms
  Grid grid;
  double nz;
  double dz = options.dz;
  double dv = options.dv;
  zlo = 0;
  zhi = options.plot_zmax;
  rlo = 0;
  rhi = options.plot_rmax;

  grid.setBounds(zlo, zhi, rhi);
  grid.setSpacing(dz, dv);
  grid.init();

  if(options.verbose) {
    cout << "Grid:" << endl;
    cout << "dz = " << dz << endl;
    cout << "dv = " << dv << endl;
    cout << "nr = " << grid.nr << endl;
    cout << "nz = " << grid.nz << endl;
    cout << "zlo = " << grid.zlo << endl;
    cout << "zhi = " << grid.zhi << endl;
    cout << "rVals = ";
    for(int i=0; i<grid.nr; i++) {
      cout << grid.rVals[i] << ", ";
    }
    cout << endl;
  }

  hDroplet = new TH2D("hDroplet","hDroplet",grid.nr,grid.rVals,grid.nz,grid.zlo,grid.zhi);
  hDroplet->SetStats(0);
  hDroplet->SetMinimum(0);
  hDroplet->SetMaximum(2.0);

  // TODO: This is probably not the way to do this.
  // Maybe pass the grid object instead.
  monolayer.hMono = new TH1D("hMono", "hMono", grid.nr, grid.rVals);

  // TODO: This is pretty messy.
  // Want to use actual coords for this one
  // Rather than shifted z coords as above
  zlo = simDataPtr->simBounds.zlo;
  zhi = simDataPtr->simBounds.zhi;
  nz = (int) ceil((zhi - zlo)/dz);
  // If dz doesn't evenly divide zhi-zlo, shift zhi up slightly.
  zhi = zlo + nz*dz;
  cout << "=======" << endl;
  cout << "Creating hLiquidDens" << endl;
  cout << "zlo = " << zlo << endl;
  cout << "zhi = " << zhi << endl;
  cout << "nz = " << nz << endl;
  cout << "dz = " << dz << endl;
  cout << "=======" << endl;
  hLiquidDens = new TH1D("hLiquidDens", "hLiquidDens", nz, zlo, zhi);
}

double Droplet::getMass() {
  // Total mass in amu
  double mass;

  // Integrate over whole cylindrical volume
  // All bins have equal volume by design
  // This quantity has units g/cm^3*AA^3
  mass = hDroplet->Integral() * options.dv;
  // Convert units to amu
  mass /= NANO_DENS_TO_MACRO;
  return mass;
}

double Droplet::getMass1D() {
  double mass;
  // TODO: Check scaling on hLiquidDens.
  // TODO: This is actually just mass in cylinder, not total
  mass = hLiquidDens->Integral() * options.dz * PI * square(rDensCyl);
  mass /= NANO_DENS_TO_MACRO;
  return mass;
}

void Droplet::dropletCalculations() {
  // TODO: Set boundary points from options
  double rBulkMax = 200.0;
  double chi2;

  // TODO: Why is monoLimits necessary here?
  // Get circle
  bulk.findBoundaryPoints(hDroplet, bulk.gCirclePoints, "a", monolayer.zlim, monolayer.hMono, rBulkMax, monolayer.radius, bulk.radius, /*unnecessary*/simDataPtr->framePtr->frameStep, /*legend*/(TLegend*)NULL, (TLine**)/*tanhLines*/NULL, (TText**)/*tanhTexts*/NULL, (TPaveText*)/*tanhTextBox*/NULL);
  chi2 = bulk.fitCircle(bulk.gCirclePoints, bulk.circle, rBulkMax, simDataPtr->framePtr->time);

  // Calculate quantities
  bulk.radius = bulk.circle.Intersect(monolayer.zlim[1]);
  bulk.contactAngle = bulk.circle.ContactAngle();
  bulk.height = bulk.circle.Height();
}
