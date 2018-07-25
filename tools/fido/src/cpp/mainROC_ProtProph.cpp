#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include "Vector.h"
#include <sstream>

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 2 )
    {
      ifstream fin(argv[1]);
      
      Array<string> protsAtThreshold;
      Array<string> cumulativeAtThreshold;
      string line;

      double prob, lastProb=-1;

      while ( fin >> prob && getline(fin, line) )
	{
	  istringstream ist(line);
	  ist >> protsAtThreshold;

	  if ( prob != lastProb && lastProb != -1.0 )
	    {
	      cout << cumulativeAtThreshold << endl;
	    }

	  cumulativeAtThreshold = concatonate( cumulativeAtThreshold, protsAtThreshold );

	  lastProb = prob;
	}
      // print the last one just in case...
      // it may be a duplicate point
      cout << cumulativeAtThreshold << endl;
    }
  else
    {
      cerr << "usage: ProtProph2ROC <xml extracted ranked list file>" << endl << endl;
    }
  return 0;
}


