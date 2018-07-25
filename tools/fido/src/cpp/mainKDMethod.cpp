#include <iostream>
#include <math.h>
#include <fstream>
#include <string>

#include "KDMethod.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 6 )
    {
      srand(0);
      cout.precision(8);

      int n = atoi(argv[1]);
      int dim = atoi(argv[2]);
      double low = atof(argv[3]);
      double high = atof(argv[4]);

      Numerical numer(atof(argv[5]));

      Array<const Vector*> points;

      int k;
      for (k=0; k<n; k++)
	{
	  Array<double> arrayRow(dim);
	  Random::fillRandomUniform(arrayRow, low, high);
	  Vector*row = new Vector(arrayRow);

	  points.add(row);
	}

      Array<Array<int> > links;
      /***
	  KDMethod slowCheck(numer, points);
      cout << "Slow: " << endl;
      links = slowCheck.slowSimilar();
      for (k=0; k<links.size(); k++)
	{
	  cout << k << ": " << links[k] << endl;
	}
      
      cout << endl;
      ***/

      KDMethod fastTest(numer, points);
      cout << endl << "Fast: " << endl;
      links = fastTest.getSimilar();
      for (k=0; k<links.size(); k++)
	{
	  links[k].sortA();
	  cout << k << ": " << links[k] << endl;
	}
      
      cout << endl;
    }
  else
    {
      cerr << "KDMethod_Test <int size> <int dim> <low> <high> <epsilon>" << endl;
    }
  return 0;
}

