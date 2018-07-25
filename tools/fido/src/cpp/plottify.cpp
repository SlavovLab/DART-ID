#include "Matrix.h"
#include <string>
#include <sstream>

int main(int argc, char**argv)
{
  if ( argc == 2 )
    {
      string both = argv[1];
      istringstream istBoth(both);

      Array<string> arrayLhs, arrayRhs;

      istBoth >> arrayLhs;
      istBoth >> arrayRhs;

      if ( arrayLhs.size() != arrayRhs.size() )
	{
	  cerr << "Error: Arrays do not have the same size" << endl;
	  exit(1);
	}

      for (int k=0; k<arrayLhs.size(); k++)
	{
	  cout << arrayLhs[k] << '\t' << arrayRhs[k] << endl;
	}
    }
  else
    {
      cerr << "usage: Plottify \"<array string> <array string>\"" << endl << endl;
    }
  return 0;
}


