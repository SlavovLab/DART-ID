#include "Matrix.h"
#include "Random.h"

#define R 1000

int main(int argc, char**argv)
{
  if ( argc == 3 )
    {
      int seed = time(NULL);
      //      seed = 0;
      //      cerr << seed << endl;
      srand(seed);

      int n = atoi(argv[1]);
      int K = atoi(argv[2]);

      int i;

      Array<string> names;

      for (i=0; i<n; i++)
	{
	  ostringstream ost;
	  ost << i;
	  names.add( "x" + ost.str() );
	}

      Array<double> xArray(n);
      xArray = Array<double>(n, 0.0001);

      Array<double> fArray(n);
      Random::fillRandomUniform(fArray, 0.0, R);


      Array<double> bArray(K);
      Random::fillRandomUniform(bArray, .01, R+.01);

      Array<Array<double> > aArray;
      for (i=0; i<K; i++)
	{
	  Array<double> aRow(n);
	  Random::fillRandomUniform(aRow, -R, R);

	  aArray.add(aRow);
	}

      cout << n << endl <<
	names << endl <<
	xArray << endl << 
	fArray << endl <<
	-bArray << endl <<
	-aArray;
    }
  else
    {
      cerr << "usage: LinProg maker <number dimensions> <number constraints>" << endl << endl;
    }
  return 0;
}


