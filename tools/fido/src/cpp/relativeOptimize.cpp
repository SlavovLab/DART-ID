#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include "Vector.h"
#include <set>


using namespace std;

int improve(const Array<double> & u, const Array<double> & v, int best)
{
  for (int k=best+1; k<u.size(); k++)
    {
      // if percent u lost is < percent v gained, accept
      double percentULost = (u[best]-u[k])/u[best];
      double percentVLost = (v[best]-v[k])/v[best];

      if ( percentULost > percentVLost )
	{
	  cout << "\ttrading when " << percentULost << " " << percentVLost << endl;
	  cout << "\timproved to @" << best << "\t" << u[best] << ", " << v[best] << endl;

	  best = k;
	}
    }

  return best;
}

int main(int argc, char**argv)
{
  if ( argc == 1 )
    {
      Array<double> u, v;
      double uEl, vEl;

      while ( cin >> uEl >> vEl )
	{
	  u.add(uEl);
	  v.add(vEl);
	}

      // you can assume that it is sorted by u (in descending order)

      // so start with the highest u element
      int best = 0;

      for (;;)
	{
	  int newBest = improve(u, v, best);
	  if ( newBest == best )
	    {
	      break;
	    }
	  best = newBest;
	  cout << "Cycle improved to @" << best << "\t" << u[best] << ", " << v[best] << endl;
	}
    }
  else
    {
      cerr << "usage: " << endl << endl;
    }
  return 0;
}


