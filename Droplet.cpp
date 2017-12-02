#include "Droplet.h"

////////////////////////
// Droplet Components //
////////////////////////

Monolayer::Monolayer()

CircularBulk::CircularBulk()
//Find the edge of droplet by tanh fitting for each row given TH2D
CircularBulk::findBoundaryPoints(TH2D* hist,TGraph *circlePointsGraph,char* aOrR,double *monoLimits,TH1D *hMono,double& rBulkMax,double &monoEdge,double &bulkEdge,int frameStep,TLegend* tanhLegend,TLine** tanhLines,TText **tanhTexts,TPaveText *tanhTextBox)
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


SphericalBulk::SphericalBulk()
CylindricalBulk::CylindricalBulk()

////////////////////
// Misc. Analysis //
////////////////////

monolayerTracker::monolayerTracker()
