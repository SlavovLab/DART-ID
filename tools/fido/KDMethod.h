#ifndef _KDMETHOD_H
#define _KDMETHOD_H

#include "Matrix.h"
#include <fstream>
#include <math.h>
#include <sstream>
#include "Random.h"

class KDMethod
{
 public:
 KDMethod(const Numerical & num, const Array<const Vector*> & vecPtrs):
  sparsify(num), orderedFeatureVectorRefs(vecPtrs), similar(vecPtrs.size())
  { 
    cout << "KDMethod linking all within " << sparsify.epsilon << endl;
    recCount = 0;
  }

  Array<Array<int> > getSimilar();
  Array<Array<int> > slowSimilar();
  
  int recCount;
 protected:

  void indent(int i) const;
  void linkSimilar(const Array<int> & from, const Array<int> & to, bool forced=false);
  void blockLinkSimilar(const Array<int> & from, const Array<int> & to);
  Array<int> blockBeginnings(const Array<double> & sortedProx);

  Vector meanVector(const Array<int> & s);

  Array<int> trueFromPoints(const Array<int> & catFromBools, int matchVal, int lowInd, int highInd);
  
  const Numerical & sparsify;
  const Array<const Vector*> & orderedFeatureVectorRefs;
  Array<Array<int> > similar;
};


#endif
