#include "Matrix.h"
#include <string>
#include <fstream>

int main(int argc, char**argv)
{
  if ( argc == 2 )
    {
      ifstream fin(argv[1]);
      Array<Array<string> > arrays;
      Array<string> array;

      for (int count = 0; count < 2; count++)
	{
	  fin >> array;
	  arrays.add(array);
	}

      int k,j;
      for (k=0; k<arrays.size()-1; k++)
	{
	  if ( arrays[k].size() != arrays[k+1].size() )
	    {
	      cerr << "Error: Arrays do not have the same size" << endl;
	      exit(1);
	    }
	}

      for (j=0; j<arrays[0].size(); j++)
	{
	  for (int k=0; k<arrays.size(); k++)
	    {
	      cout << arrays[k][j];

	      if ( k != arrays.size() - 1 )
		cout << '\t';
	    }
	  cout << endl;
	}
    }
  else
    {
      cerr << "usage: Plottify <file containing: (array string) (array string)>" << endl << endl;
    }
  return 0;
}

