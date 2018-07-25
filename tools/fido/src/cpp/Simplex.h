#ifndef _SIMPLEX_H
#define _SIMPLEX_H

#include "LPSolver.h"

class Simplex : public LPSolver
{
 public:
  Simplex(char * fname);
  Vector solve();

 protected:
  void clear();

  // matrix of constraints and variadbles
 
  void updateBasisIndices(int r, int c);
  void buildTableau();

  Matrix mat;
 
  // total number of constraints
  int K;
  // total number of constraints
  int N;

  Array<int> basisRows, basisColumns;

  void pivot(int p, int q);

};

#endif

