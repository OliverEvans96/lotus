//Read each timestep and plot 2d density from above

#include <iostream>

using namespace std;

//Count the number of water atoms in the first timestep
int countAtoms(ifstream &inFile);

//Split a string into a string array of words
string *strSplit(string str);

//Get index and coordinates from string array of words
void strToData(int *&index,double *&coords,string line,double bounds[6]);

//Read each timestep and plot 2d density from above
void topDensity(string inLoc,string outLoc)
{
	cout << "Open Stream" << endl;
	//Input Stream
	ifstream inFile;
	inFile.open(inLoc.data());
	if (inFile.good())
		cout << "Good!" << endl;
	else
		cout << "Bad!" << endl;

	cout << "Variables" << endl;
	
	//Variables
	string line;
	string input;
	int lineNum=1;
	double coords[3];
	int timestep;
	bool loopFlag1=true;
	bool loopFlag2=true;
	
	//Skip first 2 lines
	for(int i=0;i<2;i++) inFile.ignore(256,'\n');
	lineNum+=2;
	
	
	//Create Canvas
	c1 = new TCanvas();
	
	//Create Histogram
	h1 = new TH2D("h1","h1",50,-250,250,50,-250,250);
	//h1->SetStats(0);
	
	h1->SetMinimum(0.0);
	h1->SetMaximum(300.0);

	while(!inFile.eof()&&loopFlag1)
	{
		
		//Read the header
		inFile >> line >> timestep;
		inFile.ignore();
		lineNum++;
		
		cout << "Timestep " << timestep << " @ line " << lineNum-1 << endl;
		
		loopFlag2=true;
		while(loopFlag2)
		{
			getline(inFile,line);
			lineNum++;
			if(line=="") loopFlag2=false;
			else
			{
				strToData(coords,line);
				h1->Fill(coords[0],coords[1]);
			}
		}
		
		//Title
		stringstream title;
		title << "Timestep " << timestep;
		h1->SetTitle(title.str().data());
		
		//Axis Labels
		h1->GetXaxis()->SetTitle("x");
		h1->GetXaxis()->CenterTitle();
		h1->GetYaxis()->SetTitle("y");
		h1->GetYaxis()->CenterTitle();

		//Draw
		h1->Draw("colz");
		
		//Save
		stringstream filename;
		filename << outLoc << "/img/step" << setw(7) << setfill('0') << timestep << ".jpg";
		c1->SaveAs(filename.str().data());
		
		h1->Reset();
		
		//Only first timestep
 		loopFlag1=false;
	}

	delete c1;
	delete h1;
	
	inFile.close();
}

//Count the number of water atoms in the first timestep
int countAtoms(ifstream &inFile)
{
    //Count number of atoms
    bool countFlag=true;
    string line;
    int numAtoms=0;
    
    //Ignore the first 3 lines
	for(int i=0;i<3;i++) inFile.ignore(256,'\n');
    
    while(countFlag)
    {
	getline(inFile,line);
	
	//Count until reaching a line containing "TIMESTEP"
	if(line.find("TIMESTEP")!=string::npos||inFile.eof())
	{
	    countFlag=false;
		numLines-=1; //Account for the blank line
	}
	else
	    numAtoms++;
    }
	
    //Unset eof flag (if set) & return to beginning of file
    inFile.clear();
    inFile.seekg(0,ios::beg);
    
	cout << "Counted " << numAtoms << " atoms." << endl;
	
    return numAtoms;
}

//Split a string into a string array of words
string *strSplit(string str)
{
	int len=str.length();
	stringstream ss(str);
	int numWords=1;
	
	//Count number of words
	for(int ch=0;ch<len;ch++)
	if(isspace(str[ch]))
		numWords++;
	
	//Allocate array
	static string arr[16];
	
	//Read string into array
	for(int i=0;i<len;i++)
	ss >> *(arr+i);
	
	return arr;
}

//Get index and coordinates from string array of words 
void strToData(double *&coords,string line)
{
    string str;
	
	int index;
    //Split string into array of words
    string *strArr=strSplit(line);
	
    //Get index
    str=*strArr++;
    index=atoi(str.data());
    
    //Get coordinates from the second, third, and fourth elements in the array
    //string -> cstring -> double
    for(int i=0;i<3;i++)
    {
		str=*strArr++;
		*coords++=atof(str.data());
    }
}
