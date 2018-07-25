#include <iostream>
#include <math.h>
#include <fstream>
#include <string>

#include "GroupPowerBigraph.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 6 || argc == 7 )
    {
      srand(time(NULL));
      cout.precision(8);

      double gamma = atof(argv[3]);
      double alpha = atof(argv[4]);
      double beta = atof(argv[5]);

      if ( argc == 7 )
	GroupPowerBigraph::LOG_MAX_ALLOWED_CONFIGURATIONS = atof(argv[6]);

      GroupPowerBigraph gpb( argv[1], gamma, alpha, beta );
      gpb.loadKnownDetectabilitiesFast(argv[2]);

      //      gpb.showDetectabilities();
      
      //      gpb.getProteinProbs();
      //      gpb.printProteinWeights();

      //      gpb.getProteinProbs();
      //      gpb.printProteinWeights();

      //      cout << "-----------------" << endl << endl;

      gpb.getProteinProbs();
      gpb.printProteinWeights();

      //      gpb.showDetectabilities();
    }
  else
    {
      cerr << "usage: Fido <graph file> <known peptide detectabilities> <gamma> <alpha> <beta>" << endl;
      cerr << "       Fido <graph file> <known peptide detectabilities> <rbfGamma> <gamma> <alpha> <beta> <log2 of maximum number of subgraph connected states>" << endl;
    }
  return 0;
}

