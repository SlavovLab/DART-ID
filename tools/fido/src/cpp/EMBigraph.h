#ifndef _EMBigraph_h
#define _EMBigraph_h

#include "BasicGroupBigraph.h"

class EMBigraph : public BasicGroupBigraph
{
 public:

  void EM(int numberOfIterations);

 private:
  void EMIteration();
  double probR(int index);
  double expectedCountE(int index);
  double existsOtherRForE(int rIndex, int eIndex);
};

#endif
