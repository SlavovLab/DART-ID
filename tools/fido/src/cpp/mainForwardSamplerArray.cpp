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

      int seed = 0;
      //      cin >> seed;
      seed = time(NULL);
      srand(seed);

      cout << seed << endl;

      ifstream fin(argv[1]);
      ForwardSampler fs;
      fs.loadArray( fin );
      fin.close();

      fs.checkFeasibility();

      //      fs.saveArray(cout);
      //      cout << endl;

      fs.solve();
      
      fs.printResult();
    }
  else
    {
      cerr << "usage: FS-LP-Array <Dense LP-File>" << endl << endl;
    }
  return 0;
}


