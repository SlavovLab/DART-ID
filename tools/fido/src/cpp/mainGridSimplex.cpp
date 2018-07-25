#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include "GridSimplex.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 2 )
    {
      cout.precision(15);

      srand(1);

      GridSimplex gsplx(argv[1]);
      gsplx.gridSolve();
    }
  else
    {
      cerr << "usage: Grid-Simplex <LP-File>" << endl << endl;
    }
  return 0;
}


