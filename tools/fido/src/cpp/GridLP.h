#ifndef _GridLP_H
#define _GridLP_H

#include <iostream>

using namespace std;

#include "ForwardSampler.h"

class GridLP : public ForwardSampler
{
public:
 GridLP(char * fname) :
  ForwardSampler(fname)
    {
    }

  void gridSolve();

protected:
  void update(double newVal)
  {
    Array<double> bNew = ForwardSampler::b.unpack();
    bNew.back() = newVal;
    ForwardSampler::b = bNew;
  }
};

#endif

