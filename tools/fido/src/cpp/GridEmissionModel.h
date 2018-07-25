#ifndef _GridEmissionModel_H
#define _GridEmissionModel_H

#include "EmissionModel.h"
#include <iostream>

using namespace std;

class GridRange
{
 public:

  void start()
  {
    value = min;
  }

  bool inRange()
  {
    return value < max;
  }

  void advance()
  {
    value += resolution;
  }

  double getValue() const
  {
    return value;
  }

 protected:

  double value;
  double min, max;
  double resolution;
};


class GridEmissionModel : public BaseEmissionModel
{

 protected:
  GridRange grAlpha, grBeta;
  
 public:

  void start()
  {
    grAlpha.start();
  }

  
};

#endif

