#include "Parse.h"

//Masses of atoms
const double masses[3] = {26.981540,15.994915,1.007825};

int sortWaters(int argc, char* argv[])
{
	//File I/O from command line arguments
	if(argc!=3)
		cout << "2 arguments required: input & output files" << endl;

	else
	{
		//File streams
		ifstream inFile(argv[1]);
		ofstream outFile(argv[2]);

		if(!inFile.good())
		{
			cout << "Error opening " << argv[1] << endl;
			return 0;
		}

		if(!outFile.good())
		{
			cout << "Error opening " << argv[2] << endl;
			return 0;
		}

		//Count atoms
		int numAtoms=countAtoms(inFile);

		//Count timesteps
		int numSteps=countSteps(inFile);

		cout << numAtoms << " atoms and " << numSteps << " timesteps" << endl;

		//Allocate arrays
		string *header = new string[9]; //Timestep,nAtoms,Box bounds
		int *atomNumbers = new int[numAtoms];
		string *wholeLines = new string[numAtoms];

		//For each timestep
		for(int step=0;step<numSteps;step++)
		{
			//Read data into arrays
			readData(inFile,header,atomNumbers,wholeLines,numAtoms);

			//Sort arrays
			quickSort(atomNumbers,wholeLines,0,numAtoms-1);

			//Write data to output file
			writeData(outFile,header,wholeLines,numAtoms);
		}

		return 0;
	}
}

int substrateDensity(int argc, char* argv[])
{
	//File I/O from command line arguments
	if(argc!=4)
		cout << "3 arguments required: substrate input, density output, center of mass output" << endl;

	else
	{

		//File streams
		ifstream inFile(argv[1]); //Substrate input
		ofstream densFile(argv[2]); //Density output
		ofstream centerFile(argv[3]); //Center of Mass output

		if(!inFile.good())
		{
			cout << "Error opening " << argv[1] << endl;
			return 0;
		}

		if(!densFile.good())
		{
			cout << "Error opening " << argv[2] << endl;
			return 0;
		}

		if(!centerFile.good())
		{
			cout << "Error opening " << argv[3] << endl;
			return 0;
		}

		cout << "Files opened" << endl;

		//Variables
		string line;
		string arr[8];
		int timestep;
		int numAtoms=countAtoms(inFile);
		int numSteps=countSteps(inFile); //501
		double center[3];

		//z values
		//Double check this before using non-integer values for zmin,zmax,zWidth
		double zmin=-15;
		double zmax=15;
		double zWidth = 0.5;
		int nZ = (int) floor((zmax-zmin)/zWidth);

		cout << "Variables declared" << endl;

		//Allocate arrays
		double *coords = new double[3*numAtoms];
		double bounds[6];
		int *elements = new int[numAtoms];
		double *zVals = new double[nZ];
		double *density = new double[nZ];

		cout << "Arrays allocated" << endl << endl;

		//Define z-values
		for(int i=0;i<nZ;i++)
			zVals[i]=zmin+i*zWidth;

		//Write legend to first line of file
		densFile << "#TIMESTEP ";
		for(int i=0;i<nZ;i++)
			densFile << "z=" << zVals[i]	<< " ";
		densFile << endl << endl;

		//Ensure that arrays are delete []d
		try
		{
			//Loop through timesteps
			for(int step=0;step<numSteps;step++)
			{
				if(step%50==0)
					cout << "Step " << step << "/" << numSteps << endl;

				//Read data
				readStep(inFile,bounds,coords,elements,numAtoms,numSteps,timestep,arr);

				//Calculate center of mass
				centerOfMass(numAtoms,coords,elements,center);

				//Density Calculations
				calcDensity(numAtoms,bounds,coords,elements,nZ,zWidth,zVals,density,center);

				//Write data
				writeStep(densFile,centerFile,nZ,density,center,timestep);
			}
		}

		catch(...)
		{
			cout << "Catch Delete" << endl;
			//Delete arrays
			delete [] coords;
			delete [] elements;
			delete [] zVals;
			delete [] density;
			throw;
		}

		//Delete arrays
		delete [] coords;
		delete [] elements;
		delete [] zVals;
		delete [] density;

		//Close output file stream
		densFile.close();
		centerFile.close();

		cout << endl << "Calculations complete!" << endl;

		return 0;
	}
}


int calculations_centered(int argc, char* argv[])
{
	//File I/O from command line arguments
	if(argc!=3)
		cout << "2 arguments required: input & output files" << endl;

	else
	{

		//File streams
		ifstream inFile(argv[1]);
		ofstream outFile(argv[2]);

		if(!inFile.good())
		{
			cout << "Error opening " << argv[1] << endl;
			return 0;
		}

		if(!outFile.good())
		{
			cout << "Error opening " << argv[2] << endl;
			return 0;
		}

		cout << "Files opened" << endl;

		//Variables
		string line;
		int timestep;
		int numAtoms=countAtoms(inFile);
		int numMolecules=numAtoms/3;
		int numSteps=countSteps(inFile);

		cout << "Variables declared" << endl;

		//Allocate arrays
		int *OIndices = new int[numMolecules];
		double *nextOCoords = new double[3*numMolecules];
		double *OCoords = new double[3*numMolecules];
		double *prevOCoords = new double[3*numMolecules];
		int *HIndices = new int[2*numMolecules];
		double *nextHCoords = new double[2*3*numMolecules];
		double *HCoords = new double[2*3*numMolecules];
		double *velocities = new double[3*numMolecules];
		double *speeds = new double[numMolecules];
		double *dipoles = new double[3*numMolecules];

		cout << "Arrays allocated" << endl << endl;

		//Write legend to first line of file
		outFile << "Oxygen# | Coord Vector | Speed | Velocity Vector | Unit Dipole Vector" << endl << endl;

		//Ensure that arrays are delete []d
		try
		{
			//Loop through timesteps
			for(int step=0;step<=numSteps;step++)
			{
				if(step%50==0)
					cout << "Step " << step << "/" << numSteps << endl;

				//Read data
				readStep(inFile,OIndices,OCoords,prevOCoords,nextOCoords,HIndices,HCoords,nextHCoords,numMolecules,numSteps,timestep);

				//if(step<numSteps)
					//cout << "Next: (" << nextOCoords[0] << "," << nextOCoords[1] << "," << nextOCoords[2] << ")\n";


				//Don't perform any calculations for step 0
				if(step>0)
				{
					//Test
					//cout << "Current: (" << OCoords[0] << "," << OCoords[1] << "," << OCoords[2] << ")\n";

					//Centered difference for all but first and last step
					if(step<numSteps && (step!=1))
					{
						//Calculate velocities
						calculateStepVelocities(nextOCoords,prevOCoords,velocities,speeds,4e3,numMolecules);
						//cout << "Previous: (" << prevOCoords[0] << "," << prevOCoords[1] << "," << prevOCoords[2] << ")\n";
					}

					//Forward difference for first step
					else if(step==1)
					{
						//Calculate velocities
						calculateStepVelocities(nextOCoords,OCoords,velocities,speeds,2e3,numMolecules);
						//cout << "First!" << endl;
					}

					//Backward difference for last step
					else
					{
						//Calculate velocities
						calculateStepVelocities(OCoords,prevOCoords,velocities,speeds,2e3,numMolecules);
						//cout << "Previous: (" << prevOCoords[0] << "," << prevOCoords[1] << "," << prevOCoords[2] << ")\n";
						//cout << "Last!" << endl;
					}

					//cout << "Velocities: (" << velocities[0] << "," << velocities[1] << "," << velocities[2] <<  ")\n";

					//Calculate dipoles
					calculateStepDipoles(OCoords,HCoords,dipoles,numMolecules);

					//Write data
					writeStep(outFile,OIndices,OCoords,velocities,speeds,dipoles,numMolecules,timestep);
				}
			}
		}

		catch(...)
		{
			cout << "Catch Delete" << endl;
			//Delete arrays
			delete [] OIndices;
			delete [] nextOCoords;
			delete [] OCoords;
			delete [] prevOCoords;
			delete [] HIndices;
			delete [] nextHCoords;
			delete [] HCoords;
			delete [] velocities;
			delete [] speeds;
			delete [] dipoles;
			throw;
		}

		//Delete arrays
		delete [] OIndices;
		delete [] nextOCoords;
		delete [] OCoords;
		delete [] prevOCoords;
		delete [] HIndices;
		delete [] nextHCoords;
		delete [] HCoords;
		delete [] velocities;
		delete [] speeds;
		delete [] dipoles;

		//Close output file stream
		outFile.close();

		cout << endl << "Calculations complete!" << endl;

		return 0;
	}
}

//////// sortWaters

//Count the number of water atoms in the system
int countAtoms(ifstream &inFile)
{
    //Count number of atoms
    bool countFlag=true;
    string line;
    int numAtoms=0;

    //Ignore the first 9 lines
	for(int i=0;i<9;i++)
		getline(inFile,line);

    while(countFlag)
    {
		getline(inFile,line);

		if(line.find("TIMESTEP")!=string::npos||inFile.eof())
			countFlag=false;
		else
			numAtoms++;
    }

    //Unset eof flag (if set) & return to beginning of file
    inFile.clear();
    inFile.seekg(0,ios::beg);

    return numAtoms;
}

//Count the number of timesteps
int countSteps(ifstream &inFile)
{
    string line;
    int numSteps=0;
    int lineNum=0;

    //Count number of timesteps
    while(getline(inFile,line))
	{
		if(line.find("TIMESTEP")!=string::npos)
			numSteps++;
	}

    //Unset eof flag (if set) & return to beginning of file
    inFile.clear();
    inFile.seekg(0,ios::beg);

    return numSteps;
}

//Read data for one timestep from file
void readData(ifstream &inFile,string *header,int *atomNumbers,string *wholeLines,int numAtoms)
{
	bool loopFlag=true;
    string line;

	//Header
	for(int i=0;i<9;i++)
	{
		getline(inFile,line);
		*header++=line;
	}

	//Atoms section
	for(int i=0;i<numAtoms;i++)
	{
		getline(inFile,line);
		*atomNumbers++=getFirstInt(line);
		*wholeLines++=line;
	}
}

//Write sorted data to new file
void writeData(ofstream &outFile,string *header,string *wholeLines,int numAtoms)
{
	//Header
	for(int i=0;i<9;i++)
		outFile << *header++ << endl;
	//Atoms section
	for(int i=0;i<numAtoms;i++)
		outFile << *wholeLines++ << endl;
}


//Swap elements i and j of parallel arrays
void swap(int n[],string s[],int i,int j)
{
	int nt=n[i];
	string st=s[i];

	n[i]=n[j];
	s[i]=s[j];

	n[j]=nt;
	s[j]=st;
}

//Quick sort!
void quickSort(int n[],string s[],int first,int last)
{
	if(last-first>0)
	{
		//Iteration counter
		static int count=0;

		count++;
		if(count%1000000==0)
			cout << "Quick sort pass #" << count << endl;

		//Pivot
		int pivot=n[(last+first)/2];

		//Indices
		int i=first;
		int j=last;

		while(i<=j)
		{
			while(n[i]<pivot)
				i++;
			while(n[j]>pivot)
				j--;
			if(i<=j)
			{
				swap(n,s,i,j);
				i++;
				j--;
			}
		}
		quickSort(n,s,first,j);
		quickSort(n,s,i,last);
	}
}

//Get the first word from a string as an integer
int getFirstInt(string str)
{
    string firstString;
    int firstInt;

    stringstream ss(str);
    ss >> firstString;
    firstInt=atoi(firstString.data());

    return firstInt;
}

/////////// substrateDensity

//Read data for one timestep
void readStep(ifstream &inFile,double *bounds,double *coords,int *elements,int numAtoms,int numSteps,int &timestep,string *arr)
{
	//Variables
	bool loopFlag;
	string line;
	static int step=0; //Count from zero
	static int lineNum=1;

	//Read Header

	//Skip first line
	getline(inFile,line);
	//Read timestep
	inFile >> timestep;
	inFile.ignore();
	//Skip next 3 lines
	for(int i=0;i<3;i++) inFile.ignore(256,'\n');
	//Read box bounds
	for(int i=0;i<3;i++)
	{
		inFile >> bounds[2*i] >> bounds[2*i+1];
		inFile.ignore();
	}
	//Ignore next line
	inFile.ignore(256,'\n');

	lineNum+=9;

	//Read Atoms
	for(int i=0;i<numAtoms;i++)
	{
		lineNum++;
		getline(inFile,line);

		//Get next oxygen data from current line
		strToData(coords,elements,line,bounds,arr);
		//cout << "e=" << *(elements-1) << endl;
		//cout << "c=(" << *(coords-3) << "," << *(coords-4) << "," << *(coords-1) << ")" << endl;

	}

	//Increment timestep counter
	step++;

}

//Calculate density(z) for one timestep
void calcDensity(int numAtoms,double *bounds,double *coords,int *elements,int nZ,double zWidth,double *zVals,double *density,double *center)
{
	//x-y area is 171312.33624 AA^2 for spherical droplets
	double volume=(bounds[1]-bounds[0])*(bounds[3]-bounds[2])*zWidth;

	//unit conversion factor
	double convFact=10/6.022;

	int bin;

	//Start on z coordinate
	coords+=2;

	for(int i=0;i<numAtoms;i++)
	{
		//Subtract center of mass
		*coords-=center[2];

		//Identify which bin the atom falls in
		bin=(int)floor((*coords-zVals[0])/zWidth);

		//Add mass to bin
		//We are subtracting one from the element # because the file being read counts from 1, while C++ arrays count from 0.
		if(bin<nZ)
			density[bin]+=masses[*elements++-1];
		//If at the very end, add to last bin
		else
			density[bin-1]+=masses[*elements++-1];

		//Skip to next z coordinate
		coords+=3;
	}

	//Scale density from g to g/cc
	for(int i=0;i<nZ;i++)
		density[i]*=convFact/volume;
}

//Calculate the center of mass of the substrate for one timestep
void centerOfMass(int numAtoms,double *coords,int *elements,double *center)
{
	//Count total mass
	double totalMass=0;

	//Maximum x distance to consider for center of mass calculations
	double cutoff=50;

	//Reset center of mass to {0,0,0}
	for(int i=0;i<3;i++)
		center[i]=0;

	for(int j=0;j<numAtoms;j++)
	{
		//Check x position
		if(abs(*coords)<cutoff)
		{
			//Add x,y,z masses of each atom to center
			for(int i=0;i<3;i++)
				center[i]+=*coords++*masses[*elements-1];

			//Add mass of this atom to total mass
			totalMass+=masses[*elements++-1];

		}
		//Otherwise, skip to next molecule
		else
		{
			coords+=3;
			elements+=1;
		}
	}

	//Divide by total mass
	for(int i=0;i<3;i++)
		center[i]/=totalMass;
}

//Write data for one timestep
void writeStep(ofstream &densFile,ofstream &centerFile,int nZ,double *density,double *center,int timestep)
{
	//Write timestep
	densFile << timestep << " ";
	centerFile << timestep << " ";

	//Write densities
	for(int i=0;i<nZ;i++)
		densFile << *density++ << " ";

	//Write center of mass
	for(int i=0;i<3;i++)
		centerFile << *center++ << " ";

	//New line
	densFile << endl;
	centerFile << endl;
}

//Count the number of timesteps
int countSteps(ifstream &inFile,int numAtoms)
{
	cout << "Counting timesteps..." << endl;
    string line;
    int lineNum=0;
	int numSteps;

    //Count number of timesteps
    while(getline(inFile,line))
    {
		lineNum++;
		if(lineNum%100000000==0)
			cout << "line " << lineNum << endl;
    }

    //Unset eof flag (if set) & return to beginning of file
    inFile.clear();
    inFile.seekg(0,ios::beg);

    numSteps=lineNum/(numAtoms+9);
	cout << "Counted " << numSteps << " steps." << endl;
    return numSteps;
}

//Split a string into a string array of words
string *strSplit(string str,string *arr)
{
    int len=str.length();
    stringstream ss(str);
    const int numWords=8;

	/*
	//Count number of words
    for(int ch=0;ch<len;ch++)
	if(isspace(str[ch]))
	    numWords++;
	*/

    //Allocate array
    //string arr[numWords];

    //Read string into array
    for(int i=0;i<len;i++)
      ss >> *(arr+i);
}

//Get index and coordinates from string array of words
void strToData(double *&coords,int *&elements,string line,double bounds[6],string *strArr)
{
    string str;
  	string img;

    //Split string into array of words
    strSplit(line,strArr);

    //Get element
    str=*++strArr;
    *elements++=atoi(str.data());

    //Get coordinates from the third, fourth, and fifth elements in the array
    //string -> cstring -> double
    for(int i=0;i<3;i++)
    {
		str=*++strArr;
		img=*(strArr+3); //image: ix,iy,iz
		//Unscale and unwrap data
		*coords++=bounds[2*i]+(bounds[2*i+1]-bounds[2*i])*(atof(str.data())+atof(img.data()));
    }
}

///////// calculations_centered

//Read data for one timestep
void readStep(ifstream &inFile,int *OIndices,double *OCoords,double *prevOCoords,double *nextOCoords,int *HIndices,double *HCoords,double *nextHCoords,int numMolecules,int numSteps,int &timestep)
{
	//Variables
	bool loopFlag;
	string line;
	static int nextTimestep=0;
	static int step=0; //Count from zero
	static int lineNum=0;
	double bounds[6];

	//Read Header

	//Skip first line
	getline(inFile,line);
	//Save timestep #
	if(step>0)
		timestep = nextTimestep;
	inFile >> nextTimestep;
	inFile.ignore();
	//Skip next 3 lines
	for(int i=0;i<3;i++) inFile.ignore(256,'\n');
	//Read box bounds
	for(int i=0;i<3;i++)
	{
		inFile >> bounds[2*i] >> bounds[2*i+1];
		inFile.ignore();
	}
	//Ignore next line
	inFile.ignore(256,'\n');

	//Read Atoms
	for(int i=0;i<numMolecules;i++)
	{
		getline(inFile,line);
		lineNum++;

		//Previous for all but steps 0,1
		if(step>1)
		{
			//Save previous oxygen coords
			for(int i=0;i<3;i++)
				*prevOCoords++=*(OCoords+i); //previous H coords aren't important
		}

		//Current for all but step 0
		if(step>0)
		{
			//Save current oxygen and hydrogen coordinates
			//For the last step, the nextOCoords pointer is not incremented, so we've got to distinguish here.
			if (step<numSteps)
			{
				for(int k=0;k<3;k++)
				{

					*OCoords++=*(nextOCoords+k);
					for(int j=0;j<2;j++)
						*HCoords++=*(nextHCoords+2*k+j); //Need H coords for dipole calculation
				//HEY! LOOK AT THIS!
				//I just added the above j for loop. It seemed to be working without it, and I'm not seeing any difference, but it seems like it should be there. I'm too tired to understand what's goin on, though. Good luck!
				}
			}
			else
			{
				for(int k=0;k<3;k++)
				{

					*OCoords++=*nextOCoords++;
					for(int j=0;j<2;j++)
						*HCoords++=*nextHCoords++; //Need H coords for dipole calculation
				}
			}
		}

		//Next for all but last step (step numSteps)
		if(step<numSteps)
		{
			//Get next oxygen data from current line
			strToData(OIndices,nextOCoords,line,bounds);

			//Get hydrogen atom data from next two lines
			for(int j=0;j<2;j++)
			{
				//Read hydrogen data
				getline(inFile,line);
				lineNum++;
				strToData(HIndices,nextHCoords,line,bounds);
			}
		}
	}

	//Increment timestep counter
	step++;

}

//Get index and coordinates from string array of words
void strToData(int *&index,double *&coords,string line,double bounds[6])
{
  string str;
	string img;

  // TODO: This is not going to fly
  string *strArr;

  //Split string into array of words
  strSplit(line, strArr);

  //Get index
  str=*strArr++;
  *index++=atoi(str.data());

  //Get coordinates from the third, fourth, and fifth elements in the array
  //string -> cstring -> double
  for(int i=0;i<3;i++)
    {
      str=*++strArr;
      img=*(strArr+3); //image: ix,iy,iz
      *coords++=bounds[2*i]+(bounds[2*i+1]-bounds[2*i])*(atof(str.data())+atof(img.data())); //Unscale and unwrap data
    }
}

//Calculate velocities for one timestep
void calculateStepVelocities(double *nextOCoords,double *prevOCoords,double *velocities,double *speeds,double dt,int numMolecules)
{
	//Variables
	static double sumSq;

	for(int mol=0;mol<numMolecules;mol++)
	{
		//Reset sumSq
		sumSq=0;

		for(int i=0;i<3;i++)
		{
			*velocities=(*nextOCoords++ - *prevOCoords++)/2e3; //dt = 2000fs => units = AA/fs

			//Add to sumSq and move velocity pointer
			sumSq+=(*velocities * *velocities++);
		}

		//Calculate speed from sumSq
		*speeds++=sqrt(sumSq);
	}
}


//Calculate dipoles for one timestep
void calculateStepDipoles(double *OCoords,double *HCoords,double *dipoles,int numMolecules)
{
	//Coordinates of centroid between hydrogen atoms
	double avg[3];

	//non-normalized dipole
	double nonNormalized[3];

	//For magnitude
	double magnitude;
	double sumSq;

	//Loop through molecule
	for(int mol=0;mol<numMolecules;mol++)
	{
		//Reset sumSq
		sumSq=0;

		//Loop through coordinates in a molecule
		for(int i=0;i<3;i++)
		{
			//Calculate centroid
			avg[i]=(*(HCoords+i)+(*(HCoords+i+3)))/2;
			//Calculate dipole
			nonNormalized[i]=avg[i]-*OCoords++;

			//Add to sumSq
			sumSq+=nonNormalized[i]*nonNormalized[i];
		}

		//Calculate magnitude of dipole
		magnitude=sqrt(sumSq);

		//Save normalized dipole
		for(int i=0;i<3;i++)
			*dipoles++=nonNormalized[i]/magnitude;

		//Move HCoords pointer to next molecule (2 atoms, 3 coordinates each)
		HCoords+=2*3;
	}
}

//Write data for one timestep
void writeStep(ofstream &outFile,int *OIndices,double *OCoords,double *velocities,double *speeds,double *dipoles,int numMolecules,int timestep)
{
	//Mark timestep
	outFile << "TIMESTEP " << timestep << endl;

	//Loop through molecules
	for(int molecule=0;molecule<numMolecules;molecule++)
	{
		//Oxygen #
	    outFile << *OIndices++ << " ";
		//Coord Vector
	    for(int i=0;i<3;i++)
			outFile << *OCoords++ << " ";
		//Speed
		outFile << *speeds++ << " ";
		//Velocity Vector
		for(int i=0;i<3;i++)
			outFile << *velocities++ << " ";
		//Unit Dipole Vector
	    for(int i=0;i<3;i++)
			outFile << *dipoles++ << " ";
	    outFile << endl;
	}

	outFile << endl;
}
