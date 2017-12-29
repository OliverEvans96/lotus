#include "MDBase.h"

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
          numAtoms-=1; //Account for the blank line
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
