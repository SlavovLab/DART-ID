#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include "ForwardSampler.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 2 )
    {
      cout.precision(15);

      int seed;
      //      cin >> seed;
      seed = time(NULL);
      srand(seed);

      cout << seed << endl;

      //      ifstream fin(argv[1]);
      ForwardSampler fs(argv[1]);
      //      fs.load( fin );
      //      fin.close();

      fs.outputMathematica("/tmp/lp.nb");
      fs.checkFeasibility();

      //      fs.saveArray(cout);
      //      cout << endl;

      fs.solve();

      fs.printResult();
    }
  else
    {
      cerr << "usage: FS-LP <Sparse LP-File>" << endl << endl;
    }
  return 0;
}


