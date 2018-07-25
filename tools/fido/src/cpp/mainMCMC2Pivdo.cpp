#include <iostream>
#include <math.h>
#include <fstream>
#include <string>

#include "PivdoSplitter.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 3 )
    {
      cout.precision(8);

      PivdoSplitter bb;

      ifstream graph(argv[1]), pepProph(argv[2]);
      bb.readFromMCMC(graph, pepProph);

      bb.outputPivdo(cout);
    }
  else
    {
      cerr << "usage: Fido <MCMC graph file> <MCMC peptideProphet file>" << endl;
    }
  return 0;
}

