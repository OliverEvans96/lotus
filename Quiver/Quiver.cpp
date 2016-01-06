// Vector field plot from two 2D histograms

#include "Quiver.h"

// Constructors
Quiver::Quiver(int nx, double xlo, double xhi, int ny, double ylo, double yhi)
{
	// Arguments
	this->nx = nx;
	this->xlo = xlo;
	this->xhi = xhi;
	this->ny = ny;
	this->ylo = ylo;
	this->yhi = yhi;

	// Create canvas
	bool canvasCreated = false;

	// Number of data points so far
	N = 0;

	// Histograms
	h = new TH2D("h","h",nx,xlo,xhi,ny,ylo,yhi); // Dummy for plot
	hCount = new TH2I("hCount","hCount",nx,xlo,xhi,ny,ylo,yhi); // Count entries
	hvx = new TH2D("hvx","hvx",nx,xlo,xhi,ny,ylo,yhi); // v_x
	hvy = new TH2D("hvy","hvy",nx,xlo,xhi,ny,ylo,yhi); // v_y

	// Levels
	nLevels = 49;
	levels = new double[nLevels];

	// Pointers
	a = new TArrow*[nx*ny];
	l = new TText*[nx*ny];
	m = new TMarker*[nx*ny];

	// For setting color levels
	bool setLevels = false;

	// Draw Options
	drawMarkers = false;
	drawText = false;

	// Arrow Parameters
	arrowAngle = 40;
	arrowWidth = 2;
	padding = 0.05;

	// Title
	h->SetTitle("Vector Field");
}

// Destructors
Quiver::~Quiver()
{
	delete cQuiver;
	delete h;
	delete hCount;
	delete hvx;
	delete hvy;

	// Delete individual elements
	for(int i=0;i<nx*ny;i++)
	{
		delete a[i];
		delete l[i];
		delete m[i];
	}

	// Delete arrays
	delete [] levels;
	delete [] a;
	delete [] l;
	delete [] m;
}
void Quiver::Fill(double x, double y, double vx, double vy)
{
	hCount->Fill(x,y);
	hvx->Fill(x,y,vx);
	hvy->Fill(x,y,vy);

	N++;
}

void Quiver::Reset()
{

	if(canvasCreated)
	{
		cQuiver->Clear();
	}

	h->Reset();
	hCount->Reset();
	hvx->Reset();
	hvy->Reset();
	N=0;

	x1.clear();
	y1.clear();

	vxk.clear();
	vyk.clear();
	norm.clear();
	k=0;
}

void Quiver::Calculate()
{
	// Calculate norms
	k = 0;
	for (int i=0;i<nx;i++)
	{
		for (int j=0;j<ny;j++)
		{
			// Sums -> Averages (Divide by number of elements per bin)
			nElements = int(hCount->GetBinContent(i+1,j+1));

			if (nElements != 0)
			{
				vxk.push_back(hvx->GetBinContent(i+1,j+1) / nElements);
				vyk.push_back(hvy->GetBinContent(i+1,j+1) / nElements);
			}
			else
			{
				vxk.push_back(0);
				vyk.push_back(0);
			}


			x1.push_back(h->GetXaxis()->GetBinCenter(i+1));
			y1.push_back(h->GetYaxis()->GetBinCenter(j+1));

			norm.push_back(sqrt(vxk[k]*vxk[k]+vyk[k]*vyk[k]));

			k++;
		}
	}

	if(!setLevels)
	{
		minVal=minSet(norm,k);
		maxVal=maxSet(norm,k);
	}
}

void Quiver::Draw(TCanvas *c,int pad)
{
	// First, perform calculations
	Calculate();

	// Draw hist & color scale
	for (int i=1;i<nLevels;i++)
		levels[i] = minVal + (maxVal - minVal) / (nLevels - 1) * (i);
	levels[0]=1e-50;
	
	canvasCreated = true;
	cQuiver = c;
	cQuiver->cd(pad);


	// Set plot attributes
	gStyle->SetOptStat(0);
	h->SetContour(nLevels, levels);
	h->SetBinContent(0,0,0);
	h->DrawClone("col text"); // Draw axes
	h->GetZaxis()->SetRangeUser(minVal, maxVal);
	h->Draw("z same"); // Draw color scale

	// Draw arrows
	k = 0;
	for (int i=0;i<nx;i++)
	{
		for (int j=0;j<ny;j++)
		{
			minWidth = min(h->GetXaxis()->GetBinWidth(i+1),h->GetYaxis()->GetBinWidth(j+1));
			arrowLen = (1.0-padding)*minWidth;
			factor = arrowLen/(norm[k]);

			if (norm[k] != 0)
			{
				xLen = factor*vxk[k];
				yLen = factor*vyk[k];

				a[k] = new TArrow(x1[k]-xLen/2,y1[k]-yLen/2,x1[k]+xLen/2,y1[k]+yLen/2,arrowHead,"|>");
				a[k]->SetAngle(arrowAngle);
				a[k]->SetLineWidth(arrowWidth);

				// Determine color
				if(norm[k] < minVal)
					color = 51;
				else if(norm[k] > maxVal)
					color = 100;
				else
					color = 51+int((norm[k]-minVal)*nLevels/(maxVal-minVal));
				

				a[k]->SetFillColor(color);
				a[k]->SetLineColor(color);
				a[k]->Draw();

				// Convert norm to string for printing
				normStream << norm[k];
				normString = normStream.str();
				normStream.str(""); // Clear stringstream

				if(drawText)
				{
					l[k] = new TText(x1[k],y1[k],normString.data());
					l[k]->Draw();
				}

				if(drawMarkers)
				{
					m[k] = new TMarker(x1[k],y1[k],4);
					m[k]->SetMarkerStyle(4);
					m[k]->SetMarkerColor(kBlack);
					m[k]->Draw();
				}
			}
			
			k++;

		}

	}
}

void Quiver::SetLevels(double min, double max)
{
	minVal = min;
	maxVal = max;
	setLevels = true;
}

void Quiver::SetDrawOptions(bool drawMarkers,bool drawText)
{
	this->drawMarkers = drawMarkers;
	this->drawText = drawText;
}

void Quiver::SetArrowParams(double arrowAngle,double arrowHead, int arrowWidth, double padding)
{
	this->arrowAngle = arrowAngle;
	this->arrowHead = arrowHead;
	this->arrowWidth = arrowWidth;
	this->padding = padding;
}

void Quiver::SetTitle(const char *title)
{
	h->SetTitle(title);
}

void Quiver::SaveAs(const char *filename)
{
	cQuiver->SaveAs(filename);
}

double Quiver::min(double x, double y)
{
	if (x<y)
		return x;
	else
		return y;
}

double Quiver::minSet(vector<double>x,int N)
{
	double min;

	// Don't want to return 0 unless all values are 0.
	for(int i=0;i<N;i++)
	{
		if( x[i] != 0)
		{
			min = x[i];
			break;
		}
	}


	for (int i=1;i<N;i++)
	{
		if ((x[i]<min) && (x[i] != 0))
		{
			min=x[i];
		}
	}
	
	return min;
}

double Quiver::maxSet(vector<double> x,int N)
{
	double max;

	// Don't want to return 0 unless all values are 0.
	for(int i=0;i<N;i++)
	{
		if(x[i] != 0)
		{
			max = x[i];
			break;
		}
	}

	for (int i=1;i<N;i++)
	{
		if (x[i]>max)
			max=x[i];
	}
	
	return max;
}
