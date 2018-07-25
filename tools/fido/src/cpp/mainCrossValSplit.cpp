#include <iostream>
#include <math.h>
#include <fstream>
#include <string>

#include "PivdoSplitter.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 4 )
    {
      cout.precision(8);

      PivdoSplitter bb;

      ifstream fin(argv[1]);
      bb.read(fin);

      int numSplits = atoi(argv[2]);

      for (int k=0; k<numSplits; k++)
	bb.outputSplitPivdo(k, numSplits, argv[3]);
    }
  else
    {
      cerr << "usage: Fido <pivdo file> <number of splits> <dest path>" << endl;
    }
  return 0;
}

