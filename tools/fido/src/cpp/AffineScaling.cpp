#include "AffineScaling.h"

AffineScaling::AffineScaling(char * fname) :
  stopChecker(1e-10)
{
  ifstream fin(fname);

  load(fin);

  addPositivityConstraints();
  //  normalizeLP();

  A = -A;
  b = -b;
  //  f = -f;
  objective = f*x;

  gamma = .9;
}

Vector AffineScaling::solve()
{
  do
    {
      cout << "Objective " << f*x << endl;
      if ( ! singleAdvance() )
	break;
      iterations++;
    } while ( ! stoppingCriteria() );

  /***
  // removed for timing
  cout << "______________________________" << endl;
  cout << "Solved in " << cycle << " iterations" << endl;
  cout << "x is: " << x << endl;
  cout << "Final objective: " << f*x << endl;
  ***/

  return x;
}

bool AffineScaling::singleAdvance()
{
  Vector v = b - A*x;

  Vector vNegSquared( v.size() );
  for ( Set::Iterator iter = v.beginNonzero(); iter != v.endNonzero(); iter++ )
    {
      vNegSquared.addElement( *iter, 1/pow(v[ *iter ], 2.0) );
    }

  Vector hx = (A.transpose() * Matrix::diagonalMatrix(vNegSquared) * A).inverse() * f;

  Vector hv = -A*hx;

  Set st = hv < 0.0;
  if ( st.isEmpty() )
    {
      cerr << "Problem unbounded" << endl;
      exit(1);
    }

  double alpha = gamma * (-v / hv)[ st ].min();
  x = x + alpha * hx;
  lastObjective = objective;
  objective = f*x;

  return true;
}

bool AffineScaling::stoppingCriteria()
{
  //  return false;
  return stopChecker.isNonpos((objective - lastObjective) / fabs(lastObjective));
}

