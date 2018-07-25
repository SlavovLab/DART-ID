#include "GridSimplex.h"

void GridSimplex::gridSolve()
{
  double chg = -.1;
  for (double gamma = 0.0; gamma >= -15.0; gamma += chg)
    {
      //      cout << "RUNNING: " << gamma << endl << endl;
      update(gamma);

      solve();

      // performance hack for larger ones--
      // ignore after passing ROC50 interest
      if ( ( x > 0.0 ).size() > 200 )
	exit(0);

      //      cout << names[ x > 0.0 ] << endl;
    }
}
