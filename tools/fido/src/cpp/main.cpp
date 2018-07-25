#include <iostream>
#include <math.h>
#include <fstream>
#include <string>

#include "GramSchmidt.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 2 )
    {
      int N = atoi(argv[1]);
      srand(0);

      cout.precision(20);

      Array<Array<double> > matArray(N, Array<double>(N, 0.0) );
      Array<double> b(N, 0.0);

      for (int k=0; k<N; k++)
	{
	  for (int j=0; j<N; j++)
	    {
	      matArray[k][j] = rand() % 100;
	    }

	  b[k] = rand() % 10;
	}
 
      Matrix mat = matArray;

      cout << "Solution is: " << solve(mat, b).unpack() << endl;
      cout << "Fast solution is: " << fastSolve(mat, b).unpack() << endl;
   }
  else
    {
      cerr << "usage: FS-LP <Pivdo File>" << endl << endl;
    }
  return 0;
}

