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
  if ( argc == 3 )
    {
      Array<double> fpLhs, tpLhs, fpRhs, tpRhs;

      ifstream finLhs(argv[1]);
      ifstream finRhs(argv[2]);

      finLhs >> fpLhs >> tpLhs;
      finRhs >> fpRhs >> tpRhs;

      Array<double> fpInter, tpDiffInter;
      for (double t=0.0; t<=50; t+=1)
	{
	  fpInter.add(t);
	  tpDiffInter.add( ROCInterpolator::interpolate(fpLhs, tpLhs, t) - ROCInterpolator::interpolate(fpRhs, tpRhs, t) );
	}

      cout << fpInter << " " << tpDiffInter << endl;
    }
  else
    {
      cerr << "usage: ROCDifferenceInterpolate <ROC file> <ROC file>" << endl;
    }
  return 0;
}


