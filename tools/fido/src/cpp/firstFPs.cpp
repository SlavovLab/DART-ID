#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include "Vector.h"

using namespace std;

#define STOP_THRESH 3

int matchCount( Array<string> positiveNames, Array<string> cumulativeAtThreshold )
{
  int count = 0;
  
  for (int k=0; k<cumulativeAtThreshold.size(); k++)
    {
      for (int j=0; j<positiveNames.size(); j++)
	{
	  if ( cumulativeAtThreshold[k] == positiveNames[j] )
	    count++;
	}
    }

  return count;
}

void printFPs( Array<string> falsePositiveNames, Array<string> cumulativeAtThreshold )
{
  Array<string> fpAtThreshold;
  for (int k=0; k<cumulativeAtThreshold.size(); k++)
    {
      for (int j=0; j<falsePositiveNames.size(); j++)
	{
	  if ( cumulativeAtThreshold[k] == falsePositiveNames[j] )
	    fpAtThreshold.add( cumulativeAtThreshold[k] );
	}
    }

  cout << fpAtThreshold << endl;
}

int main(int argc, char**argv)
{
  if ( argc == 2 )
    {
      ifstream fin(argv[1]);
      Array<string> truePositiveNames, falsePositiveNames;
      fin >> truePositiveNames;
      fin >> falsePositiveNames;

      Array< Array<double> > tpAtFp(5000);

      Array<string> cumulativeAtThreshold;
      string line;

      Array<double> fpPoints, tpPoints;


      while ( getline(cin, line) )
	{
	  istringstream ist(line);
	  ist >> cumulativeAtThreshold;

	  //	  cout << "Read: " << cumulativeAtThreshold << endl;

	  //	  int tp = matchCount(truePositiveNames, cumulativeAtThreshold);
	  int fp = matchCount(falsePositiveNames, cumulativeAtThreshold);

	  if ( fp >= STOP_THRESH )
	    {
	      printFPs(falsePositiveNames, cumulativeAtThreshold);
	      return 0;
	    }
	}

    }
  else
    {
      cerr << "usage: First-FPs <file with TPs FPs>  <  string lists file" << endl << endl;
    }
  return 0;
}


