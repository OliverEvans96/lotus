//////////////////////////////////////////////////////////////////////////
//            Line Integral Convolution for Flow Visualization          //
//                            Initial  Version                          //
//                              May 15, 1999                            //
//                                   by                                 //
//                            Zhanping Liu                              //
//                     (zhanping@erc.msstate.edu)                       //
//                               while with                             //
//                          Graphics  Laboratory                        //
//                           Peking  University                         //
//                              P. R.  China                            //
//----------------------------------------------------------------------//
//                             Later  Condensed                         //
//                             May 4,  2002                             //
//          VAIL (Visualization Analysis & Imaging Laboratory)          //
//                  ERC  (Engineering Research Center)                  //
//                     Mississippi State University                     //
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
//This code was developed  based on  the original algorithm  proposed by//
//Brian Cabral  and  Leith (Casey) Leedom  in the paper  "Imaging Vector//
//Fields Using Line Integral Convolution", published  in  Proceedings of//
//ACM  SigGraph 93,  Aug 2-6,  Anaheim,  California,  pp. 263-270, 1993.//
//Permission to use, copy, modify, distribute and sell this code for any//
//purpose is hereby granted without fee,  provided that the above notice//
//appears in all copies  and  that both that notice and  this permission//
//appear in supporting documentation. The  developer  of this code makes//
//no  representations  about  the  suitability  of  this  code  for  any//
//purpose.  It is provided "as is"  without express or implied warranty.//
//////////////////////////////////////////////////////////////////////////

// Refer to Enhanced LIC:
// Enhanced Line Integral Convolution with Flow Feature Detection (1997)
// by Arthur Okada , David Lane
// In SPIE Vol. 3017 Visual Data Exploration and Analysis IV
// http://citeseerx.ist.psu.edu/viewdoc/download;jsessionid=62D238FE3BCA070245D7748174523482?doi=10.1.1.35.6799&rep=rep1&type=pdf

#include "FieldViz.h"

// Constructor
FieldViz::FieldViz()
{
    licLength = 40;
    nIter = 2;
    histEqExp = 0.85;

    strcpy(cmapFile, "phase_cmap.txt");

    DISCRETE_FILTER_SIZE = 2048;
    LINE_SQUARE_CLIP_MAX = 100000.0;
    VECTOR_COMPONENT_MIN = 0.050000 ;

    p_LUT0 = (double*) malloc( sizeof(double) * DISCRETE_FILTER_SIZE);
    p_LUT1 = (double*) malloc( sizeof(double) * DISCRETE_FILTER_SIZE);
}

FieldViz::FieldViz(int n_xPix, int n_yPix)
{
    licLength = 40;
    nIter = 2;
    histEqExp = 0.85;

    strcpy(cmapFile, "phase_cmap.txt");

    DISCRETE_FILTER_SIZE = 2048;
    LINE_SQUARE_CLIP_MAX = 100000.0;
    VECTOR_COMPONENT_MIN = 0.050000;

    p_LUT0 = (double*) malloc( sizeof(double) * DISCRETE_FILTER_SIZE);
    p_LUT1 = (double*) malloc( sizeof(double) * DISCRETE_FILTER_SIZE);

	SetNumPixels(n_xPix,n_yPix);
	Init();
}

FieldViz::FieldViz(int n_xPix, int n_yPix, double xmin, double xmax, double ymin, double ymax)
{
    licLength = 40;
    nIter = 2;
    histEqExp = 0.85;

    strcpy(cmapFile, "phase_cmap.txt");

    DISCRETE_FILTER_SIZE = 2048;
    LINE_SQUARE_CLIP_MAX = 100000.0;
    VECTOR_COMPONENT_MIN = 0.050000;

    p_LUT0 = (double*) malloc( sizeof(double) * DISCRETE_FILTER_SIZE);
    p_LUT1 = (double*) malloc( sizeof(double) * DISCRETE_FILTER_SIZE);

	SetNumPixels(n_xPix,n_yPix);
	Init();
	SetBounds(xmin,xmax,ymin,ymax);
}

// Destructor
FieldViz::~FieldViz()
{
    free(pVectr);    pVectr = NULL;
    free(p_LUT0);    p_LUT0 = NULL;
    free(p_LUT1);    p_LUT1 = NULL;
    free(pNoise);    pNoise = NULL;
    free(pImage);    pImage = NULL;
    free(cmap);    cmap = NULL;
    free(vecMags); vecMags = NULL;
    free(vecAngles); vecMags = NULL;
    free(colorIndices); colorIndices = NULL;


    // Delete field arrays
    for(int ii=0; ii<n_xPix; ii++)
    {
        delete [] xField[ii];
        delete [] yField[ii];
        delete [] xFieldTmp[ii];
        delete [] yFieldTmp[ii];
    }
    delete [] xField;
    delete [] yField;
    delete [] xFieldTmp;
    delete [] yFieldTmp;
}

void FieldViz::SetColormap(char* _cmapFile)
{
      strcpy(cmapFile, _cmapFile);
}

void FieldViz::SetNumPixels(int _n_xPix, int _n_yPix)
{
    n_xPix = _n_xPix;
    n_yPix = _n_yPix;
}

// Initialize after setting number of pixels
void FieldViz::Init()
{
    pNoise = (unsigned char* ) malloc( sizeof(unsigned char) * n_xPix * n_yPix     );
    pImage = (unsigned char* ) malloc( sizeof(unsigned char) * n_xPix * n_yPix     );
    pVectr = (double* ) malloc( sizeof(double) * n_xPix * n_yPix * 2);
    cmap = (double* ) malloc( sizeof(double) * 256 * 3 );
    vecMags = (double* ) malloc( sizeof(double) * n_xPix * n_yPix );
    vecAngles = (double* ) malloc( sizeof(double) * n_xPix * n_yPix );
    colorIndices = (int*) malloc( sizeof(int) * n_xPix * n_yPix);

    // Allocate 2D field arrays
    xField = new double*[n_xPix];
    yField = new double*[n_xPix];
    xFieldTmp = new double*[n_xPix];
    yFieldTmp = new double*[n_xPix];
    for(int ii=0; ii<n_xPix; ii++)
    {
        xField[ii] = new double[n_yPix];
        yField[ii] = new double[n_yPix];
        xFieldTmp[ii] = new double[n_yPix];
        yFieldTmp[ii] = new double[n_yPix];
    }

    LoadColormap(cmapFile);
    MakeWhiteNoise(pNoise);
    GenBoxFiltrLUT(DISCRETE_FILTER_SIZE);
}

// Prepare vector field for LIC
void FieldViz::PreProcess()
{
	DivideFieldByFrames();
    NormalizVectrs();
	NormField();
    //CalculateColor(vecMags);
	CalculateColor(vecAngles);
}

// Set number of iterations and exponent for histogram equalization
void FieldViz::SetLICParams(int _nIter, double _histEqExp)
{
    nIter = _nIter;
    histEqExp = _histEqExp;
}

void FieldViz::PerformLIC()
{
    // Repeat LIC as many times as desired
    for(int iter=0; iter<nIter; iter++)
	{
		cout << "LIC iteration " << iter+1 << endl;
        FlowImagingLIC(pNoise, pImage, licLength);
	}
}

// Improve output image quality after LIC
void FieldViz::PostProcess()
{
    HighPass();
    HistEqual(histEqExp);
}

// Provide 2d atoms & velocities to be turned into an evenly sampled vec. field.
void FieldViz::SetData(double* _xx, double* _yy, double* _vx, double* _vy, int _nAtoms)
{
    xx = _xx;
    yy = _yy;
    vx = _vx;
    vy = _vy;
    nAtoms = _nAtoms;

	cout << "Set FieldViz data with " << nAtoms << " atoms" << endl;
}

// Set all field values
void FieldViz::ResetField()
{
	int index;

	// Reset number of frames
	nFrames++;

	// Reset grid & arrays
	for(int ii=0; ii<n_xPix; ii++)
	{
		for(int jj=0; jj<n_yPix; jj++)
		{
			xField[ii][jj] = 0;
			yField[ii][jj] = 0;

			index = ii * n_yPix + jj;
			pVectr[ii] = 0;
		}
	}

	// Reset number of frames 
	nFrames = 0;
}

// Specify the number of atoms
void FieldViz::SetNAtoms(int _nAtoms)
{
	nAtoms = _nAtoms;
}

// Provide 3d Cartesian atom positions & velocities
// Convert to cylindrical coordinates, ignoring azimuthal dependence
// xx, yy, zz - positions
// vx, vy, vz - velocities
// nAtoms - number of atoms
void FieldViz::SetCylDataFromCart(double* _xx, double* _yy, double* zz,
                                  double* _vx, double* _vy, double* vz,
                                  int nAtoms)
{
	// Radial position & velocity
    double rr[nAtoms], vr[nAtoms];

    // Loop over atoms & convert to cylindrical coordinates
    for(int ii=0; ii<nAtoms; ii++)
    {
        rr[ii] = sqrt(_xx[ii]*_xx[ii] + _yy[ii]*_yy[ii]);

        // Calculate velocity w/ chain rule
        // Remember: x = r*cos(th) => cos(th) = x/r. Same for sin(th).
        vr[ii] = ( _vx[ii]*_xx[ii] + _vy[ii]*_yy[ii] ) / rr[ii];
    }

    // Save 2d cylindrical data, ignoring theta
    SetData(rr, zz, vr, vz, nAtoms);
}


// Set bounds for visualization
void FieldViz::SetBounds(double _xMin, double _xMax, double _yMin, double _yMax)
{
    xMin = _xMin;
    xMax = _xMax;
    yMin = _yMin;
    yMax = _yMax;

    xRange = xMax - xMin;
    yRange = yMax - yMin;

    // Pixel resolution - data units per pixel
    xPixRes = xRange / n_xPix;
    yPixRes = yRange / n_yPix;
}

// Load colormap from file
void FieldViz::LoadColormap(char* _cmapFile)
{
    strcpy(cmapFile, _cmapFile);
    ifstream cmapStream(cmapFile);
    int index;

    if(!cmapStream.good())
        cout << "Failed to load cmap @ '" << cmapFile << "'" << endl;

    for(int i=0; i<256; i++)
    {
        index = 3*i;
        cmapStream >> cmap[index] >> cmap[index+1] >> cmap[index+2];
    }

    cmapStream.close();
}

// Calculate color index from magnitude
void FieldViz::CalculateColor(double* vecMags)
{
    int index;
    double mag;
    double min = MinVal(vecMags, n_xPix*n_yPix);
    double max = MaxVal(vecMags, n_xPix*n_yPix);
    double range = max - min;

    for(int i=0; i< n_xPix; i++)
    {
        for(int j=0; j<n_yPix; j++)
        {
            index = j*n_xPix + i;
            mag = vecMags[index];
            if(mag > min && mag < max)
                colorIndices[index] = (int) floor(256*(mag - min)/range);
            else if (mag <= min)
                colorIndices[index] = 0;
            else
                colorIndices[index] = 255;
        }
    }
}

// Find minimum value from array
double FieldViz::MinVal(double* array, int n)
{
    double min = array[0];

    for(int i=0; i<n; i++)
        if(array[i] < min)
            min = array[i];

    return min;
}

// Find minimum value from array
double FieldViz::MaxVal(double* array, int n)
{
    double max = array[0];

    for(int i=0; i<n; i++)
        if(array[i] > max)
            max = array[i];

    return max;
}

// Write vector field to file in 2 files: xField, yField
// Each file contains a 2D array over space
// First file contains x components, 2nd has y components
// FILES ARE WRITTEN IN MATRIX FORMAT, NOT XY FORMAT
// format = "ij" : write in matrix form, matrix element (0,0) in resulting file corresponding to (0,N-1) in cartesian grid
// format = "xy" : write in cartesian form, matrix element (0,0) in resulting file corresponding to (0,0) in cartesian grid
void FieldViz::WriteField(const char* xFile, const char* yFile, const char* format)
{
	ofstream xStream(xFile);
	ofstream yStream(yFile);

	int ww = 15;
	int pp = 5;

	int index, i_prime, j_prime;

	xStream.precision(pp);
	yStream.precision(pp);
	xStream << scientific;
	yStream << scientific;

	if(strcmp(format,"ij")==0)
	{
		for(int ii=0; ii<n_yPix; ii++)
		{
			for(int jj=0; jj<n_xPix; jj++)
			{
				i_prime = jj;
				j_prime = n_yPix - 1 - ii;

				xStream << setw(ww) << xField[i_prime][j_prime];
				yStream << setw(ww) << yField[i_prime][j_prime];
			}
			xStream << endl;
			yStream << endl;
		}
	}

	else if(strcmp(format,"xy")==0)
	{
		for(int ii=0; ii<n_xPix; ii++)
		{
			for(int jj=0; jj<n_yPix; jj++)
			{
					xStream << setw(ww) << xField[ii][jj];
					yStream << setw(ww) << yField[ii][jj];
			}
			xStream << endl;
			yStream << endl;
		}
	}

	else
		cout << "Invalid format for WriteField" << endl;

	xStream.close();
	yStream.close();
}

// Write output to one 4-column file:
// x, y, vx, vy
void FieldViz::WriteField4Col(char* outFile)
{
	ofstream outStream(outFile);

	int ww = 15;
	int pp = 5;

	double xx, yy;

	outStream.precision(pp);
	outStream << scientific;

	for(int ii=0; ii<n_xPix; ii++)
	{
		for(int jj=0; jj<n_yPix; jj++)
		{
			xx = xMin + xPixRes * ii;
			yy = yMin + yPixRes * jj;
			outStream 
				<< setw(ww) << xx
				<< setw(ww) << yy
				<< setw(ww) << xField[ii][jj]
				<< setw(ww) << yField[ii][jj]
				<< endl;
		}
	}

	outStream.close();
}

// Add new frame to previous field (Needs to be divided later)
void FieldViz::AddFieldFromData(double eps)
{
	// Copy data to temporary field arrays
	for(int ii=0; ii<n_xPix; ii++)
	{
		for(int jj=0; jj<n_yPix; jj++)
		{
			xFieldTmp[ii][jj] = xField[ii][jj];
			yFieldTmp[ii][jj] = yField[ii][jj];
		}
	}
	
	// Generate field with new data
	GenFieldFromData(eps);

	// Add previous field to current
	for(int ii=0; ii<n_xPix; ii++)
	{
		for(int jj=0; jj<n_yPix; jj++)
		{
			xField[ii][jj] += xFieldTmp[ii][jj];
			yField[ii][jj] += yFieldTmp[ii][jj];
		}
	}

	// Increment number of frames since last reset
	nFrames++;
}

// Divide field by number of frames included in field
void FieldViz::DivideFieldByFrames()
{
	for(int ii=0; ii<n_xPix; ii++)
	{
		for(int jj=0; jj<n_yPix; jj++)
		{
			xField[ii][jj] /= nFrames;
			yField[ii][jj] /= nFrames;
		}
	}
}

// Generate evenly sampled vector field from position & velocity data
// eps - radius of ball over which to average velocities for each Pixel
void FieldViz::GenFieldFromData(double eps)
{
    // Index of block containing pixel
    int xBlockInd, yBlockInd;
	// Index of block presently being searched for nearby atoms
	int xSearchBlock, ySearchBlock;

	// Dummy indices
	int index, i_prime, j_prime;

    // Pixel location in data coordinates
    double xPixLoc, yPixLoc;

    // Quadrant of block: +1 for upper, -1 for lower
    // in each direction
    int xQuad, yQuad;

    // Pointer to index list for current block
    vector<int>* p_curBlockInd;

    // Index of atom relative to overall position/velocity vectors
    int atomIndex;
    // Atom data coordinates
    double xAtom, yAtom;

    // Number of atoms in neighborhood of pixel
    int n_nearbyAtoms;

	// Diameter of circles of influence
	double diam = 2*eps;

	// Distance between atom & pixel & resulting weight
	double dist, weight;

    ///////////////////////////////////////////////////////////
    // Divide spatial domain into square blocks of size diam //
    ///////////////////////////////////////////////////////////

    // Determine number of blocks in each direction
    n_xBlocks = (int) ceil(xRange / diam);
    n_yBlocks = (int) ceil(yRange / diam);

    // Allocate vector of atom indices for each block
    indVecs = new vector<int>*[n_xBlocks];
    for(int ii=0; ii<n_xBlocks; ii++)
        indVecs[ii] = new vector<int>[n_yBlocks];

    // Loop through atoms & place in appropriate blocks
    for(int kk=0; kk<nAtoms; kk++)
    {
        xBlockInd = (int) floor(xx[kk] / diam);
        yBlockInd = (int) floor(yy[kk] / diam);
		//cout << "xBlockInd: " << xBlockInd << " / " << n_xBlocks << endl;
		//cout << "yBlockInd: " << yBlockInd << " / " << n_yBlocks << endl;

        // Save indices of atoms to bin vectors
		// Ignore atoms outside of domain
		if(0 <= xBlockInd && xBlockInd < n_xBlocks &&
		   0 <= yBlockInd && yBlockInd < n_yBlocks)
		{
			indVecs[xBlockInd][yBlockInd].push_back(kk);
		}
    }

    ///////////////////////////////////////////////////////////////
    // Loop through pixels & calculate neighborhood avg velocity //
    ///////////////////////////////////////////////////////////////

    // Loop through pixels
    for(int ii=0; ii<n_xPix; ii++)
    {
        for(int jj=0; jj<n_yPix; jj++)
        {
            // Reset field components for this pixel
            xField[ii][jj] = 0;
            yField[ii][jj] = 0;

            // Calculate location of center of pixel in data coordinates
            xPixLoc = (ii + .5) * xPixRes;
            yPixLoc = (jj + .5) * yPixRes;

            // Determine which block this pixel belongs to
            xBlockInd = (int) floor(xPixLoc / diam);
            yBlockInd = (int) floor(yPixLoc / diam);

            // Consider the present block to be broken up into quadrants
            // Which quadrant the pixel falls in determines the adjacent
            // blocks which need to be searched to find all atoms within
            // a ball of radius eps from the pixel.

            // xQuad / yQuad are +1 for upper quad., -1 for lower quad.
            // Outer floor is solely there to supress compiler warning
            xQuad = (int) floor( 2 * floor( 2.0/diam * fmod(xPixLoc,diam) ) - 1 );
            yQuad = (int) floor( 2 * floor( 2.0/diam * fmod(yPixLoc,diam) ) - 1 );

            // Reset counter of atoms near pixel
            n_nearbyAtoms = 0;

            // Search through appropriate adjacent blocks for nearby atoms
            for(int xBool=0; xBool<2; xBool++)
            {
                for(int yBool=0; yBool<2; yBool++)
                {

					// Block we're presently searching through
					xSearchBlock = xBlockInd + xBool*xQuad;
					ySearchBlock = yBlockInd + yBool*yQuad;

					// Only consider this block if it's inside of the domain
					// Could happen if pixel is near boundary
					if( 0 <= xSearchBlock && xSearchBlock < n_xBlocks &&
						0 <= ySearchBlock && ySearchBlock < n_yBlocks)
					{
						// Point to index list for presently considered block
						p_curBlockInd = &indVecs[xSearchBlock][ySearchBlock];

						// Loop through atoms in block
						for(int kk=0; kk<p_curBlockInd->size(); kk++)
						{
							// Get overall index of current atom
							atomIndex = (*p_curBlockInd)[kk];

							// Extract current atom data coordinates
							xAtom = xx[atomIndex];
							yAtom = yy[atomIndex];

							// Check distance from pixel center
							dist = sqrt(  (xPixLoc - xAtom)*(xPixLoc - xAtom)
								    	+ (yPixLoc - yAtom)*(yPixLoc - yAtom));
							if(dist < eps)
							{
								// Add velocity to average
								// Weighted by distance from pixel
								weight = SmoothStep(dist,eps);
								xField[ii][jj] += vx[atomIndex]*weight;
								yField[ii][jj] += vy[atomIndex]*weight;

								// Add atom to counter
								n_nearbyAtoms++;
							}
						}
					}
				}
            }

            // Divide by number of nearby atoms found (if any)
			if(n_nearbyAtoms > 0)
			{
				xField[ii][jj] /= n_nearbyAtoms;
				yField[ii][jj] /= n_nearbyAtoms;
			}

			// Rotate & Transpose from XY to IJ
			i_prime = n_yPix - jj - 1;
			j_prime = ii;
			index = n_xPix * i_prime + j_prime;

            // Save field to single pointer array
            //pVectr[2*(ii*n_yPix+jj)] = xField[ii][jj];
            //pVectr[2*(ii*n_yPix+jj)+1] = yField[ii][jj];
			pVectr[2*index] = xField[ii][jj];
			pVectr[2*index+1] = xField[ii][jj];


        }
    }

    // Delete index vector array
    for(int ii=0; ii<n_xBlocks; ii++)
        delete [] indVecs[ii];
    delete [] indVecs;
}

// 3x3 High pass filter
void FieldViz::HighPass()
{
    // Sharpening convolution kernel
    double K[9] = {    0,   -1./5,    0,
                   -1./5,       5,-1./5,
                       0,   -1./5,    0};

    double* newimg = (double*) malloc(sizeof(double) * n_xPix * n_yPix);
    int outer_ind, inner_ind, ind1, ind2;
    double tmpdbl;

    // Loop through image pixels
    for(int i=0; i<n_xPix; i++)
    {
        for(int j=0; j<n_yPix; j++)
        {
            outer_ind = j*n_xPix+i;
            newimg[outer_ind] = 0;
            // Loop through 3x3 convolution kernel
            // & adjacent pixels in image
            for(int m=0; m<3; m++)
            {
                for(int n=0; n<3; n++)
                {
                    // Extend boundary
                    ind1 = i-1+m;
                    ind2 = j-1+n;
                    if(ind1 < 0)
                        ind1 = 0;
                    else if (ind1 > n_xPix-1)
                        ind1 = n_xPix-1;
                    if(ind2 < 0)
                        ind2 = 0;
                    else if (ind2 > n_yPix-1)
                        ind2 = n_yPix-1;
                    inner_ind = ind2*n_xPix + ind1;
                    tmpdbl = (double) pImage[inner_ind];
                    newimg[outer_ind] += tmpdbl * K[3*n+m];
                }
            }
        }
    }

    // Rescale to [0,256]
    double min = MinVal(newimg, n_xPix*n_yPix);
    double max = MaxVal(newimg, n_xPix*n_yPix);
    double range = max-min;
    for(int i=0; i<n_xPix*n_yPix; i++)
        newimg[i] = 256 * (newimg[i]-min)/range;

    // Update image
    for(int i=0; i<n_xPix*n_yPix; i++)
        pImage[i] = (unsigned char) floor(newimg[i]);

    free(newimg);
}

// Approximate histogram equalization filter by raising to power
void FieldViz::HistEqual(double pwr)
{
    cout << "Histogram equalization: pwr=" << pwr << endl;
    int ind;
    double tmp;

    for(int i=0; i<n_xPix; i++)
    {
        for(int j=0; j<n_yPix; j++)
        {
            ind = n_xPix*j + i;
            tmp = (double) pImage[ind];
            // Transform to [-1,1]
            tmp = 2 * (tmp/256. - .5);
            // Raise to power to push closer to -1 or 1
            if(tmp >= 0)
                tmp = pow(tmp,pwr);
            else
                tmp = -pow(fabs(tmp),pwr);
            // Convert back to [0,256]
            tmp = 256 * (tmp/2 + .5);
            pImage[ind] = (unsigned char) floor(tmp);
        }
    }
}

// Weight function for distance in GenFieldFromData
// Smoothstep-like function
// S(1)=2/a
// S(a)=0
// S'(0)=0
// S'(a)=0
// int(S,{0,a}) = 1
double FieldViz::SmoothStep(double tt, double aa)
{
	return 2/aa * (2*(tt/aa)*(tt/aa)*(tt/aa) - 3*(tt/aa)*(tt/aa) + 1);
}

///        synthesize a saddle-shaped vector field     ///
void FieldViz::SyntheszSaddle()
{
    for(int  j = 0;  j < n_yPix;  j ++)
    for(int  i = 0;  i < n_xPix;  i ++)
    {
        int     index = (  (n_yPix - 1 - j) * n_xPix + i  )  <<  1;
        pVectr[index    ] = - ( j / (n_yPix - 1.0f) - 0.5f );
        pVectr[index + 1] =     i / (n_xPix - 1.0f) - 0.5f;
    }
}


///        normalize the vector field     ///
void FieldViz::NormalizVectrs()
{
        for(int     j = 0;  j < n_yPix;  j ++)
        for(int     i = 0;  i < n_xPix;  i ++)
        {
            int        index = (j * n_xPix + i) << 1;
            double    vcMag = double(  sqrt( double(pVectr[index] * pVectr[index] + pVectr[index + 1] * pVectr[index + 1]) )  );
            // Save vector magnitude to array
            vecMags[index >> 1] = vcMag;

            double    scale = (vcMag == 0.0f) ? 0.0f : 1.0f / vcMag;
            pVectr[index    ] *= scale;
            pVectr[index + 1] *= scale;
        }
}

void FieldViz::NormField()
{
	double mag;
	int index, i_prime, j_prime;

	for(int ii=0; ii<n_xPix; ii++)
	{
		for(int jj=0; jj<n_yPix; jj++)
		{
			mag = sqrt(xField[ii][jj]*xField[ii][jj]
					 + yField[ii][jj]*yField[ii][jj]);

			// Save magnitude
			j_prime = ii;
			i_prime = n_yPix - 1 - jj;

			index = i_prime*n_xPix + j_prime;
			vecMags[index] = mag;

			// Save angle
			vecAngles[index] = atan2(xField[ii][jj],yField[ii][jj]);

			// Normalize vector
			if (mag > 0)
			{
				xField[ii][jj] /= mag;
				yField[ii][jj] /= mag;
			}
		}
	}
}
				


///        make white noise as the LIC input texture     ///
void FieldViz::MakeWhiteNoise(unsigned char* _pNoise)
{
    for(int  j = 0;   j < n_yPix;  j ++)
    for(int  i = 0;   i < n_xPix;  i ++)
    {
        int  r = rand();
        r = (  (r & 0xff) + ( (r & 0xff00) >> 8 )  ) & 0xff;
        _pNoise[j * n_xPix + i] = (unsigned char) r;
    }
}


///        generate box filter LUTs     ///
void FieldViz::GenBoxFiltrLUT(int LUTsiz)
{
   for(int  i = 0;  i < LUTsiz;  i ++) p_LUT0[i] = p_LUT1[i] = i;
}


///        write the LIC image to a PPM file     ///
void FieldViz::WriteImage2PPM(const char* f_name)
{
    cout << "Write to file: '" << f_name << "'" << endl;
    double dblval;
    unsigned char color[3];
    int index;
    FILE*    o_file;
    if(   ( o_file = fopen(f_name, "w") )  ==  NULL   )
    {
        printf("Can't open output file\n");
        return;
    }

    fprintf(o_file, "P6\n%d %d\n255\n", n_xPix, n_yPix);

    for(int  j = 0;  j < n_yPix;  j ++)
    for(int  i = 0;  i < n_xPix;  i ++)
    {
        index = j*n_xPix + i;
        dblval = pImage[index];
        for(int k=0; k<3; k++)
            color[k] = (unsigned char) floor(dblval * cmap[3*colorIndices[index] + k]);

        fprintf(o_file, "%c%c%c", color[0], color[1], color[2]);
    }

    fclose (o_file);    o_file = NULL;
}


///        flow imaging (visualization) through Line Integral Convolution     ///
void FieldViz::FlowImagingLIC(unsigned char*  pNoise, unsigned char*  pImage,double krnlen )
{
        int        vec_id;                        ///ID in the VECtor buffer (for the input flow field)
        int        advDir;                        ///ADVection DIRection (0: positive;  1: negative)
        int        advcts;                        ///number of ADVeCTion stepS per direction (a step counter)
        int        ADVCTS = int(krnlen * 3);    ///MAXIMUM number of advection steps per direction to break dead loops
		int		   xindex, yindex;

        double    vctr_x;                        ///x-component  of the VeCToR at the forefront point
        double    vctr_y;                        ///y-component  of the VeCToR at the forefront point
        double    clp0_x;                        ///x-coordinate of CLiP point 0 (current)
        double    clp0_y;                        ///y-coordinate of CLiP point 0    (current)
        double    clp1_x;                        ///x-coordinate of CLiP point 1 (next   )
        double    clp1_y;                        ///y-coordinate of CLiP point 1 (next   )
        double    samp_x;                        ///x-coordinate of the SAMPle in the current pixel
        double    samp_y;                        ///y-coordinate of the SAMPle in the current pixel
        double    tmpLen;                        ///TeMPorary LENgth of a trial clipped-segment
        double    segLen;                        ///SEGment   LENgth
        double    curLen;                        ///CURrent   LENgth of the streamline
        double    prvLen;                        ///PReVious  LENgth of the streamline
        double    W_ACUM;                        ///ACcuMulated Weight from the seed to the current streamline forefront
        double    texVal;                        ///TEXture VALue
        double    smpWgt;                        ///WeiGhT of the current SaMPle
        double    t_acum[2];                    ///two ACcUMulated composite Textures for the two directions, perspectively
        double    w_acum[2];                    ///two ACcUMulated Weighting values   for the two directions, perspectively
        double*    wgtLUT = NULL;                ///WeiGhT Look Up Table pointing to the target filter LUT
          double    len2ID = (DISCRETE_FILTER_SIZE - 1) / krnlen;    ///map a curve LENgth TO an ID in the LUT

        ///for each pixel in the 2D output LIC image///
        for(int  j = 0;     j < n_yPix;  j ++)
        for(int  i = 0;     i < n_xPix;  i ++)
        {
            ///init the composite texture accumulators and the weight accumulators///
            t_acum[0] = t_acum[1] = w_acum[0] = w_acum[1] = 0.0f;

            ///for either advection direction///
            for(advDir = 0;  advDir < 2;  advDir ++)
            {
                ///init the step counter, curve-length measurer, and streamline seed///
                advcts = 0;
                curLen = 0.0f;
                clp0_x = i + 0.5f;
                clp0_y = j + 0.5f;

                ///access the target filter LUT///
                wgtLUT = (advDir == 0) ? p_LUT0 : p_LUT1;

                ///until the streamline is advected long enough or a tightly  spiralling center / focus is encountered///
                while( curLen < krnlen && advcts < ADVCTS )
                 {
                    ///access the vector at the sample///
                    vec_id = ( int(clp0_y) * n_xPix + int(clp0_x) ) << 1;
                    vctr_x = pVectr[vec_id    ];
                    vctr_y = pVectr[vec_id + 1];

					xindex = int(clp0_x);
					yindex = n_yPix - 1 - int(clp0_y);
					vctr_x = xField[xindex][yindex];
					vctr_y = yField[xindex][yindex];

                    ///in case of a critical point///
                    if( vctr_x == 0.0f && vctr_y == 0.0f )
                    {
                        t_acum[advDir] = (advcts == 0) ? 0.0f : t_acum[advDir];           ///this line is indeed unnecessary
                        w_acum[advDir] = (advcts == 0) ? 1.0f : w_acum[advDir];
                        break;
                    }

                    ///negate the vector for the backward-advection case///
                    vctr_x = (advDir == 0) ? vctr_x : -vctr_x;
                    vctr_y = (advDir == 0) ? vctr_y : -vctr_y;

                    ///clip the segment against the pixel boundaries --- find the shorter from the two clipped segments///
                    ///replace  all  if-statements  whenever  possible  as  they  might  affect the computational speed///
                    segLen = LINE_SQUARE_CLIP_MAX;
                    segLen = (vctr_x < -VECTOR_COMPONENT_MIN) ? ( int(     clp0_x         ) - clp0_x ) / vctr_x : segLen;
                    segLen = (vctr_x >  VECTOR_COMPONENT_MIN) ? ( int( int(clp0_x) + 1.5f ) - clp0_x ) / vctr_x : segLen;
                    segLen = (vctr_y < -VECTOR_COMPONENT_MIN) ?
                             (      (    (  tmpLen = ( int(     clp0_y)          - clp0_y ) / vctr_y  )  <  segLen    ) ? tmpLen : segLen      )
                            : segLen;
                    segLen = (vctr_y >  VECTOR_COMPONENT_MIN) ?
                             (      (    (  tmpLen = ( int( int(clp0_y) + 1.5f ) - clp0_y ) / vctr_y  )  <  segLen    ) ? tmpLen : segLen      )
                            : segLen;

                    ///update the curve-length measurers///
                    prvLen = curLen;
                    curLen+= segLen;
                    segLen+= 0.0004f;
					
					//cout << "segLen = " << segLen << endl;

                    ///check if the filter has reached either end///
                    segLen = (curLen > krnlen) ? ( (curLen = krnlen) - prvLen ) : segLen;

                    ///obtain the next clip point///
                    clp1_x = clp0_x + vctr_x * segLen;
                    clp1_y = clp0_y + vctr_y * segLen;

                    ///obtain the middle point of the segment as the texture-contributing sample///
                    samp_x = (clp0_x + clp1_x) * 0.5f;
                    samp_y = (clp0_y + clp1_y) * 0.5f;

                    ///obtain the texture value of the sample///
                    texVal = pNoise[ int(samp_y) * n_xPix + int(samp_x) ];

                    ///update the accumulated weight and the accumulated composite texture (texture x weight)///
                    W_ACUM = wgtLUT[ int(curLen * len2ID) ];
                    smpWgt = W_ACUM - w_acum[advDir];
                    w_acum[advDir]  = W_ACUM;
                    t_acum[advDir] += texVal * smpWgt;

                    ///update the step counter and the "current" clip point///
                    advcts ++;
                    clp0_x = clp1_x;
                    clp0_y = clp1_y;

					//cout << "len = " << curLen << "/" << krnlen << endl;
                    ///check if the streamline has gone beyond the flow field///
                    if( clp0_x < 0.0f || clp0_x >= n_xPix || clp0_y < 0.0f || clp0_y >= n_yPix)  break;
                }
             }

            ///normalize the accumulated composite texture///
            texVal = (t_acum[0] + t_acum[1]) / (w_acum[0] + w_acum[1]);

			// Don't alter the image if no advection occured
			//if(advcts==0)
			//	texVal = pNoise[j*n_xPix + i];
			//if(advcts!=0)
			//	cout << "ADV: " << advcts << endl;

			//cout << "(" << i << "," << j << "): a=" << advcts << endl;

            ///clamp the texture value against the displayable intensity range [0, 255]
            texVal = (texVal <   0.0f) ?   0.0f : texVal;
            texVal = (texVal > 255.0f) ? 255.0f : texVal;
            pImage[j * n_xPix + i] = (unsigned char) texVal;
        }
}

