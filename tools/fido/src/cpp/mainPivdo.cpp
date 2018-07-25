#include <iostream>
#include <math.h>
#include <fstream>
#include <string>

#include "PivdoBigraph.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 2 )
    {
      cout.precision(20);

      PivdoBigraph pb;
      ifstream fin(argv[1]);
      fin >> pb;

      pb.getProteinWeights();
      pb.printProteinWeights();
    }
  else
    {
      cerr << "usage: FS-LP <Pivdo File>" << endl << endl;
    }
  return 0;
}

