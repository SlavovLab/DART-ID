#include <iostream>
#include <math.h>
#include <fstream>
#include <string>

#include "GroupPowerBigraph.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 8 || argc == 9 )
    {
      srand(time(NULL));
      cout.precision(8);

      double rbfGamma = atof(argv[3]);
      double gamma = atof(argv[4]);
      double alpha = atof(argv[5]);
      double beta = atof(argv[6]);

      int iters = atoi(argv[7]);

      if ( argc == 9 )
	GroupPowerBigraph::LOG_MAX_ALLOWED_CONFIGURATIONS = atof(argv[8]);

      GroupPowerBigraph gpb( argv[1], argv[2], rbfGamma, gamma, alpha, beta );
      
      //      gpb.getProteinProbs();
      //      gpb.printProteinWeights();

      //      gpb.getProteinProbs();
      //      gpb.printProteinWeights();

      //      cout << "-----------------" << endl << endl;

      gpb.getProteinProbs();
      //      gpb.printProteinWeights();

      gpb.printPeptideDetectabilities();
      //      cout << "-----------------" << endl;

      for (int k=0; k<iters; k++)
	{
	  cerr << "Starting iteration " << k << endl;
	  gpb.refreshModels();

	  
	  cerr << "Models have been refreshed" << endl;
	  gpb.getProteinProbs();
	  cerr << "Got protein probs" << endl;

	  //	  gpb.printProteinWeights();

	  //	  if ( k != iters - 1 )
	  //	    cerr << endl << endl << "-----------------" << endl << endl;
	}
      //      gpb.showDetectabilities();

      //      gpb.printProteinWeights();
      //      gpb.refreshModels();
      
      //      cout << "recalc. prot probs." << endl;

    }
  else
    {
      cerr << "usage: Fido <graph file> <peptide detectability file> <rbfGamma> <gamma> <alpha> <beta> <number iterations>" << endl;
      cerr << "       Fido <graph file> <peptide detectability file> <rbfGamma> <gamma> <alpha> <beta> <number iterations> <log2 of maximum number of subgraph connected states>" << endl;
    }
  return 0;
}

