#ifndef _ROCInterpolator_H
#define _ROCInterpolator_H

#include <iostream>
#include "Array.h"

using namespace std;

class ROCInterpolator
{
 public:
  static double interpolate(const Array<double> & fps, const Array<double> & tps, double fpVal);
  static double interpolateCollection(const Array<Array<double> > & fpCollection, const Array<Array<double> > & tpCollection, double fpVal);
};

#endif

