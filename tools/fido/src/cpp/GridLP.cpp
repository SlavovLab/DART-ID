#include "GridLP.h"

void GridLP::gridSolve()
{
  double chg = -.1;
  for (double gamma = 0.0; gamma >= -15.0; gamma += chg)
    {
      //      cout << "RUNNING: " << gamma << endl << endl;
      update(gamma);

      solve();
      Set st = x > 0.0;
      cout << names[ st ] << endl;

      // note: hack: bail out when there are too many points to matter for
      // the ROC50
      if ( st.size() > 100 )
	exit(0);
    }
}
