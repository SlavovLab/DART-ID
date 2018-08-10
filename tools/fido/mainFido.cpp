#include <iostream>
#include <math.h>
#include <fstream>
#include <string>

#include "GroupPowerBigraph.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 5 || argc == 6 )
    {
      srand(time(NULL));
      cout.precision(8);

      double gamma = atof(argv[2]);
      double alpha = atof(argv[3]);
      double beta = atof(argv[4]);

      if ( argc == 6 )
	GroupPowerBigraph::LOG_MAX_ALLOWED_CONFIGURATIONS = atof(argv[5]);

      cerr << "\tLoading..." << endl;
      GroupPowerBigraph gpb( argv[1], gamma, alpha, beta);

      cerr << "\tMarginalizing..." << endl;
      //      cout << gpb.logLikelihoodAlphaBetaGivenD() << endl;

      //      cout << gpb.getLogNumberStates() << endl;
      //      exit(0);

      gpb.getProteinProbs();

      //      gpb.printGraphs();
      //      cout << endl << "----------------------" << endl << endl;

      gpb.printProteinWeights();

      //      gpb.showDetectabilities();

      //      gpb.printPeptideDetectabilities();
    }
  else
    {
      cerr << "usage: Fido <graph file> <gamma> <alpha> <beta>" << endl;
      cerr << "       Fido <graph file> <gamma> <alpha> <beta> <log2 of maximum number of subgraph connected states>" << endl;
    }
  return 0;
}

