#ifndef _AffineScaling_H
#define _AffineScaling_H

#include "LPSolver.h"

using namespace std;


class AffineScaling : public LPSolver
{
public:
  AffineScaling(char * fname);
  Vector solve();

private:
  bool singleAdvance();
  bool stoppingCriteria();

  double objective, lastObjective;

  double gamma;

  Numerical stopChecker;
};

#endif

