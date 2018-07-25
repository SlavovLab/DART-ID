#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include "GridLP.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 2 )
    {
      cout.precision(15);

      srand(1);

      GridLP glp(argv[1]);
      glp.gridSolve();
    }
  else
    {
      cerr << "usage: Grid-LP <LP-File>" << endl << endl;
    }
  return 0;
}


