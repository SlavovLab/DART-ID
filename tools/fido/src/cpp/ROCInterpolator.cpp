#include "ROCInterpolator.h"

double ROCInterpolator::interpolate(const Array<double> & fps, const Array<double> & tps, double fpVal)
{
  bool success = false;

  int k;
  for (k=0; k<fps.size(); k++)
    {
      if ( fps[k] > fpVal )
	{
	  success = true;
	  break;
	}
    }

  if ( ! success )
    {
      cerr << "Problem: no place to interpolate here (value = " << fpVal << " )" << endl;
      exit(1);
    }

  // k is the first index that passes fpVal
  if ( k!= 0 )
    return tps[k-1] + (tps[k] - tps[k-1])/(fps[k]-fps[k-1])*(fpVal-fps[k-1]);
  else
    return tps[0];
}

double ROCInterpolator::interpolateCollection(const Array<Array<double> > & fpCollection, const Array<Array<double> > & tpCollection, double fpVal)
{
  double tot = 0.0;

  for (int k=0; k<fpCollection.size(); k++)
    {
      tot += interpolate(fpCollection[k], tpCollection[k], fpVal);
    }

  return tot / fpCollection.size();
}

