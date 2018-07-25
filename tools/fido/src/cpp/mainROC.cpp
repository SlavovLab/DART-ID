#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include "Vector.h"

using namespace std;

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

	  int tp = matchCount(truePositiveNames, cumulativeAtThreshold);
	  int fp = matchCount(falsePositiveNames, cumulativeAtThreshold);

	  cerr << fp << " " << tp << endl;

	  tpAtFp[ fp ].add((double) tp);

	  if ( fp > 50 )
	    break;
	}

      if ( tpAtFp[0].size() == 0 )
	{
	  fpPoints.add(0.0);
	  tpPoints.add(0.0);
	}

      int best = 0;
      for (int k=0; k<1000; k++)
	{
	  best = best > Vector(tpAtFp[ k ]).max() ? best : (int)Vector(tpAtFp[ k ]).max();

	  if ( tpAtFp[k].size() > 0 )
	    {
	      //	      cout << k << " " << Vector( tpAtFp[ k ] ).max() << endl;
	      fpPoints.add(k);
	      tpPoints.add( Vector( tpAtFp[ k ] ).max() );
	    }
	  //	  else
	  //	    cout << k << " " << best << endl;
	}

      int overallTP = truePositiveNames.size();
      int overallFP = falsePositiveNames.size();
      if ( fpPoints.back() < overallFP )
	{
	  fpPoints.add(overallFP);
	  tpPoints.add(overallTP);
	}
      else
	{
	  //	  cerr << "Declining to add big point" << endl;
	}

      cout << fpPoints << " " << tpPoints << endl;
    }
  else
    {
      cerr << "usage: <True Positives and Overall Array File>" << endl << endl;
    }
  return 0;
}


