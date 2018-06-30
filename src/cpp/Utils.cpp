#include "Utils.h"

//Count number of lines in file
int countLines(ifstream &inFile) {
  int pos = inFile.tellg();
  string line;
  int numLines=0;

  //Beginning of file
  inFile.seekg(0,ios::beg);

  while(inFile.good()) {
    getline(inFile,line);
    numLines++;
  }

  inFile.clear();
  inFile.seekg(pos);

  return numLines;
}

//Elementwise vector addition
vector<double> add(vector<double> u,vector<double> v) {
  vector<double> w=u;
  int n=v.size();
  for(int i=0;i<n;i++)
    w[i]+=v[i];
  return w;
}

//Elementwise vector multiplication
vector<double> mult(vector<double> u,vector<double> v) {
  vector<double> w=u;
  int n=v.size();
  for(int i=0;i<n;i++)
    w[i]*=v[i];
  return w;
}

double sum(vector<double> v) {
  double sum=0;
  int n=v.size();
  for(int i=0;i<n;i++)
    sum+=v[i];
  return sum;
}

//Calculate standard deviation
double stddev(vector<double> x) {
  double m = mean(x);
  double sum=0;
  for(int i=0;i<x.size();i++)
    sum+=(m-x[i])*(m-x[i]);
  return sqrt(sum/(x.size()-1));
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

double min(double *v, int n) {
  double min=v[0];

  for(int i=0;i<n;i++)
    {
      if(v[i]<min)
        min=v[i];
    }

  return min;
}

double max(double *v, int n) {
  double max=v[0];

  for(int i=0;i<n;i++)
    {
      if(v[i]>max)
        max=v[i];
    }

  return max;
}

double mean(double *v, int n) {
  int m = 0;
  for(int i=0; i<n; i++) {
    m += v[i];
  }
  return m/n;
}

double stddev(double *v, int n) {
  double m = mean(v, n);
  double sum=0;
  for(int i=0;i<n;i++)
    sum+=(m-v[i])*(m-v[i]);
  return sqrt(sum/(n-1));
}

//Square
double square(double x) {
  return x*x;
}

//Arctanh
double atanh(double x) {
  return log((1+x)/(1-x))/2;
}

bool isLess(int a,int b) {
  return (a<b);
}

//Check whether x is in v
bool isIn(int x, vector<int>v) {
  return binary_search(v.begin(),v.end(),x,isLess);
}

// Does the string contain a substring?
bool isIn(string substr, string str) {
  return str.find(substr)!=string::npos;
}

//Check whether a file exists and is not a directory
bool file_exists(const char* pathname) {
  struct stat sb;
  return (stat(pathname, &sb) == 0 && !S_ISDIR(sb.st_mode));
}

// Check whether a directory exists
// https://stackoverflow.com/questions/3828192/checking-if-a-directory-exists-in-unix-system-call
bool dir_exists(const char* pathname) {
  struct stat sb;
  return (stat(pathname, &sb) == 0 && S_ISDIR(sb.st_mode));
}

//round up to nearest multiple of a number
int roundUp(int numToRound, int multiple) {
  if (multiple == 0)
    return numToRound;

  int remainder = numToRound % multiple;
  if (remainder == 0)
    return numToRound;

  return numToRound + multiple - remainder;
}


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

//Given a TH1D and a two bin numbers, draw a line between the points and solve for where y=yc (y-0.5)
double solveLinear(TH1D *hist,int bin1,int bin2,double yc)
{
  //Find x & y values at bins
  double x1=hist->GetBinCenter(bin1);
  double x2=hist->GetBinCenter(bin2);
  double y1=hist->GetBinContent(bin1);
  double y2=hist->GetBinContent(bin2);

  //Equation of line
  double m = (y2-y1)/(x2-x1);
  double b = y1-m*x1;

  //Solve
  double xc = (yc-b)/m;
  return xc;
}

// Path operations

void stripTrailingSlash(char* strippedPath, const char* path) {
  int len = strlen(path);

  strcpy(strippedPath, path);
  if(strippedPath[len-1] == '/') {
    // Replace with null terminator
    strippedPath[len-1] = 0;
  }
}

void joinPath(char* path, const char* prefix, const char* suffix) {
  char strippedPath[256];
  // Remove trailing slash from prefix if present
  stripTrailingSlash(strippedPath, prefix);
  sprintf(path, "%s/%s", strippedPath, suffix);
}
