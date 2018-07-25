#ifndef _LPBigraph_H
#define _LPBigraph_H

#include "BasicBigraph.h"

using namespace std;

class LPBigraph : public BasicBigraph
{
public:
  LPBigraph();

  void outputLP() const;
  void sumProteins();

protected:
};

#endif

