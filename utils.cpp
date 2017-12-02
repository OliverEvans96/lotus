#include "utils.h"

//Split a string into a string vector of words
vector<string> strSplit(string str)
{
  int len=str.length();
  stringstream ss(str);
  int numWords=1;
  bool foundSpace=false;

  //Count number of words
  for(int ch=0;ch<len;ch++)
    {
      if(isspace(str[ch]))
        foundSpace=true;
      if(!isspace(str[ch])&&foundSpace)
        {
          numWords++;
          foundSpace=false;
        }
    }

  //Allocate array
  vector<string> arr(numWords);

  //Read string into array
  for(int i=0;i<len;i++)
    ss >> arr[i];

  return arr;
}

void strToData(double *coords,double *velocities,double &dipole,string line)
{
  string str;

  int index;
  //Split string into array of words
  vector<string> strArr=strSplit(line);

  //Get index
  str=strArr[0];
  index=atoi(str.data());

  //Save values
  //string -> cstring -> double
  for(int i=0;i<3;i++)
    {
      //Get coordinates from the second, third, and fourth elements in the array
      str=strArr[1+i];
      *coords++=atof(str.data());

      //Get velocities from the sixth, seventh, and eighth elements in the array
      str=strArr[5+i];
      *velocities++=atof(str.data());
    }

  //Save dipole moment cos(theta)
  dipole=atof(strArr[10].data());
}

//Choose the highest value in a vector
double max(vector<double> v)
{
  int n=v.size();
  double max=v[0];

  for(int i=0;i<n;i++)
    {
      if(v[i]>max)
        max=v[i];
    }

  return max;

}

//Find the mean of a double vector
double mean(vector<double> v)
{
  double m=0;
  for(vector<double>::iterator it = v.begin(); it != v.end(); ++it)
    m+=*it;
  m/=v.size();

  return m;
}

//Square
inline double square(double x) {return x*x;}

//Arctanh
inline double atanh(double x) {return log((1+x)/(1-x))/2;}

bool isLess(int a,int b) { return (a<b); }

//Check whether x is in v
bool isIn(int x, vector<int>v)
{
  return binary_search(v.begin(),v.end(),x,isLess);
}

//Check whether a file exists
bool file_exists(const string& name);

//round up to nearest multiple of a number
int roundUp(int numToRound, int multiple);

double findMinimum(vector<double> v)
{
  double min=v[0];
  for(int i=0;i<v.size();i++)
    {
      if(v[i]<min)
        {
          min=v[i];
        }
    }
  return min;
}

double findMaximum(vector<double> v)
{
  double max=v[0];
  for(int i=0;i<v.size();i++)
    {
      if(v[i]>max)
        {
          max=v[i];
        }
    }
  return max;
}
