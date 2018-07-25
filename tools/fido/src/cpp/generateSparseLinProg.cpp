#include "Matrix.h"
#include "Random.h"

#define R 1000

int main(int argc, char**argv)
{
  if ( argc == 4 )
    {
      int seed = time(NULL);
      //      seed = 0;
      //      cerr << seed << endl;
      srand(seed);

      int n = atoi(argv[1]);
      int K = atoi(argv[2]);
      double probabilityOfZeros = atof( argv[3] );

      int i;

      Array<string> names;

      for (i=0; i<n; i++)
	{
	  ostringstream ost;
	  ost << i;
	  names.add( "x" + ost.str() );
	}

      Array<double> xArray(n);
      xArray = Array<double>(n, 1e-5);

      Array<double> fArray(n);
      Random::fillRandomUniform(fArray, 0, R);


      Array<double> bArray(K);
      Random::fillRandomUniform(bArray, 10, R+10);

      Array<Array<double> > aArray;
      for (i=0; i<K; i++)
	{
	  Array<double> aRow(n);
	  Random::fillRandomUniform(aRow, -R, R);

	  // fill with zeros
	  for (int j=0; j<n; j++)
	    {
	      if ( Random::uniform( 0.0, 1.0 ) < probabilityOfZeros )
		{
		  aRow[j] = 0.0;
		}
	    }

	  aArray.add(aRow);
	}

      cout << n << endl <<
	names << endl <<
	Vector(xArray) << endl << 
	Vector(fArray) << endl <<
	-Vector(bArray) << endl <<
	-Matrix(aArray);
    }
  else
    {
      cerr << "usage: SparseLinProg maker <number dimensions> <number constraints> <probability of zeros>" << endl << endl;
    }
  return 0;
}


