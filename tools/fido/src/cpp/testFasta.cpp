#include <iostream>
#include <math.h>
#include <fstream>
#include <string>

using namespace std;

#include "FastaReader.h"
#include "Proteomics.h"

int main(int argc, char**argv)
{
  if ( argc == 2 )
    {
      FastaReader fr(argv[1]);

      //      Proteomics::setMinMaxLength(argv[1], 1.0);

      for ( int k=0; k<fr.names.size(); k++ )
	{
	  Array<string> peptides = Proteomics::trypsinDigest( fr.sequences[k] );

	  for (int j=0; j<peptides.size(); j++)
	    {
	      //	      if ( Proteomics::observablePeptide( peptides[j] ) )
		{
		  cout << "e " << peptides[j] << endl
		       << "r " << fr.names[k] << endl
		       << "p 0.0" << endl;
		}
	    }
	}
    }
  else
    {
      cerr << "usage: FastaTest <fasta file>" << endl << endl;
    }
  return 0;
}

