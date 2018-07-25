#include <iostream>
#include <math.h>
#include <fstream>
#include <string>

#include "GroupPowerBigraph.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 9 || argc == 5 )
    {
      srand(time(NULL));
      cout.precision(8);

      double alpha_min, alpha_max, beta_min, beta_max, resolution, gamma;

      if ( argc == 9 )
	{
	  //	  cout << "Using user defined alpha range, beta range" << endl;
	  resolution = atof( argv[3] );
	  gamma = atof( argv[4] );
	  alpha_min = atof( argv[5] );
	  alpha_max = atof( argv[6] );
	  beta_min = atof( argv[7] );
	  beta_max = atof( argv[8] ); 
	}
      else
	{
	  //	  cout << "Using user defined resolution only, std [0,1]x[0,1] alpha and beta range" << endl;
	  resolution = atof( argv[3] );
	  gamma = atof( argv[4] );

	  alpha_min = beta_min = resolution;
	  alpha_max = beta_max = 1.0;
	}

//       cout << "alpha inRange " << alpha_min << " " << alpha_max << endl;
//       cout << "beta inRange " << beta_min << " " << beta_max << endl;
//       cout << "res " << resolution << endl << endl;
//       cout << "gamma " << gamma << endl;

      GroupPowerBigraph gpb( RealRange(alpha_min, resolution, alpha_max), RealRange(beta_min, resolution, beta_max), gamma );

      ifstream graph(argv[1]), pepProph(argv[2]);
      gpb.readFromMCMC(graph, pepProph);

      cerr << "Log of total number of states: " << gpb.getLogNumberStates() << endl;

      //      gpb.getProteinProbsOverAllAlphaBeta();

      //      gpb.printProteinWeights();

    }
  else
    {
      cerr << "usage: Fido <MCMC graph file> <MCMC peptideProphet file> <resolution> <gamma>" << endl;
      cerr << "usage: Fido <MCMC graph file> <MCMC peptideProphet file> <resolution> <gamma> <alpha_min> <alpha_max> <beta_min> <beta_max>" << endl << endl;
    }
  return 0;
}

