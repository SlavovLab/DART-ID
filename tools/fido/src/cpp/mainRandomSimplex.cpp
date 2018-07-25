#include <iostream>
#include <math.h>
#include <fstream>
#include <string>

#include "Simplex.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 2 )
    {
      srand(time(NULL));
      cout.precision(20);

      Simplex simp(argv[1]);
      simp.randomSolve();
    }
  else
    {
      cerr << "usage: Simplex <Array LP File>" << endl << endl;
    }
  return 0;
}

