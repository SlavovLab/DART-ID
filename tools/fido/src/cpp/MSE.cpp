#include "Matrix.h"
#include "Vector.h"
#include <string>
#include <sstream>

const double POWER = 4;

double normBetweenPoints(double x1, double y1, double x2, double y2)
{
  double slope = (y2-y1)/(x2-x1);
  
  
}

int main(int argc, char**argv)
{
  if ( argc == 3 )
    {
      double THRESH = atof(argv[1]);

      string both = argv[2];
      istringstream istBoth(both);

      Array<double> arrayLhs, arrayRhs;

      istBoth >> arrayLhs;
      istBoth >> arrayRhs;

      if ( arrayLhs.size() != arrayRhs.size() )
	{
	  cerr << "Error: Arrays do not have the same size" << endl;
	  exit(1);
	}

      Vector v = Vector(arrayLhs) - Vector(arrayRhs);

      double tot = 0.0;
      int count = 0;
      for (int k=0; k<v.size()-1; k++)
	{
	  tot += normBetweenPoints;
	  if ( arrayLhs[k] >= THRESH )
	    {
	      break;
	    }
	}

      double greatest = min( 0.1, arrayLhs[k] );

      cout << tot / ( greatest - arrayLhs[0] );
      //      cout << pow( tot / v.size(), 1/POWER);
    }
  else
    {
      cerr << "usage: MSE FDR_threshold \"<array double> <array double>\"" << endl << endl;
    }
  return 0;
}


