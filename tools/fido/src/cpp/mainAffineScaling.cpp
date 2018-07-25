#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include "AffineScaling.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 2 )
    {
      cout.precision(15);

      int seed = 0;
      //      cin >> seed;

      //      seed = time(NULL);
      srand(seed);

      cout << seed << endl;

      AffineScaling primal(argv[1]);
      primal.solve();
      
      primal.printResult();
    }
  else
    {
      cerr << "usage: AS-LP-Array <Dense LP-File>" << endl << endl;
    }
  return 0;
}


