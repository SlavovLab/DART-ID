#include <iostream>
#include <math.h>
#include <fstream>
#include <string>

#include "EMBigraph.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 2 )
    {
      srand(time(NULL));
      cout.precision(6);

      EMBigraph emb;
      ifstream fin(argv[1]);
      emb.read(fin);

      emb.EM(100);
      emb.printProteinWeights();
    }
  else
    {
      cerr << "usage: EM <Pivdo File>" << endl;
    }
  return 0;
}

