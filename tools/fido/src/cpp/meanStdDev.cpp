#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include "Vector.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 2 )
    {
      Array<double> array;

      istringstream ist(argv[1]);

      ist >> array;
      Vector vec = array;

      cout << vec.average() << "\t" << sqrt(vec.variance()) << endl;
    }
  else
    {
      cerr << "usage:" << endl << endl;
    }
  return 0;
}


