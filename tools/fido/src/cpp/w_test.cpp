#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include "Vector.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 3 )
    {
      Array<double> lhs, rhs;

      istringstream istL(argv[1]);
      istringstream istR(argv[2]);

      istL >> lhs;
      istR >> rhs;

      //      cout << "\t" << lhs << endl;
      //      cout << "\t" << rhs << endl;

      //      cout << (Vector(lhs) - Vector(rhs)).unpack();

      Vector diff = Vector(lhs) - Vector(rhs);
      Array<double> arrayDiff = diff.unpack();

      vector<pair<double, int> > sorter(arrayDiff.size());

      int k;
      for (k=0; k<arrayDiff.size(); k++)
	{
	  sorter[k] = pair<double,int>( fabs(arrayDiff[k]), k);
	  //	  cout << sorter[k].second << " " << sorter[k].first << endl;
	}

      sort(sorter.begin(), sorter.end());
      
      //      cout << "Sorted" << endl;
      //      for (k=0; k<arrayDiff.size(); k++)
      //	{
      //	  cout << sorter[k].second << " " << sorter[k].first << endl;
      //	}      

      double W = 0.0;

      // note: this does not yet resolve ties properly...

      for (k=0; k<int(sorter.size()); k++)
	{
	  int ind = sorter[k].second;

	  //	  if ( arrayDiff[ind] > 0.0 )
	  //	    W += (k+1);
	  if ( arrayDiff[ind] < 0.0 )
	  	    W -= (k+1);
	}

      cout << "n = " << sorter.size() << endl;
      cout << "W = " << W << endl;

      //      cout << diff.unpack() << endl;
      //      cout << diff.average() << " " << sqrt(diff.variance()) << endl;

      //      cout << diff.average() / sqrt(diff.variance()) << endl;
    }
  else
    {
      cerr << "usage:" << endl << endl;
    }
  return 0;
}


