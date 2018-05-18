#include "Droplet.h"

////////////////////////
// Droplet Components //
////////////////////////

Monolayer::Monolayer() {};

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
void Monolayer::findMonoLimits(TH1D *hWaterDens,double *monoLimits)
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

CircularBulk::CircularBulk() {};

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
double CircularBulk::fitCircle(TGraph* g,CircleFit &circle,double xMax,int timestep)
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

SphericalBulk::SphericalBulk() {};

CylindricalBulk::CylindricalBulk() {};


////////////////////
// Misc. Analysis //
////////////////////
