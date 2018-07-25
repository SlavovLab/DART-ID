#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include "Vector.h"
#include "ROCInterpolator.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 2 )
    {
      ifstream fin(argv[1]);
      Array<Array<double> > fpCollection, tpCollection;
      
      string buff;
      while ( getline(fin, buff) )
	{
	  istringstream ist(buff);
	  Array<double> fp, tp;
	  ist >> fp >> tp;

	  fpCollection.add(fp);
	  tpCollection.add(tp);
	}

      Array<double> fpInter, tpInter;
      for (double t=0.0; t<=50; t+=1)
	{
	  fpInter.add(t);
	  tpInter.add( ROCInterpolator::interpolateCollection(fpCollection, tpCollection, t) );
	}

      cout << fpInter << " " << tpInter << endl;
    }
  else
    {
      cerr << "usage: <True Positives and Overall Array File>" << endl << endl;
    }
  return 0;
}


