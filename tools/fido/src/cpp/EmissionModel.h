#ifndef _EmissionModel_H
#define _EmissionModel_H

#include <math.h>
#include "Combinatorics.h"

class EmissionModel
{
 protected:
  double alpha, beta;
 public:
  double associatedEmission() const;
  double spontaneousEmission() const;
  double probabilityNoEmissionFrom(int numActiveProts) const;
};

#endif

