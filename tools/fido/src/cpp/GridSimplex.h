#ifndef _GridSimplex_H
#define _GridSimplex_H

#include <iostream>

using namespace std;

#include "Simplex.h"

class GridSimplex : public Simplex
{
public:
 GridSimplex(char * fname) :
  Simplex(fname)
    {
    }

  void gridSolve();

protected:
  void update(double newVal)
  {
    Array<double> bNew = b.unpack();
    bNew.back() = newVal;
    b = bNew;
    
    buildTableau();
  }
};

#endif

