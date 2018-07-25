#include "PivdoBigraph.h"

void PivdoBigraph::getProteinWeights()
{
  for (int k=0; k<proteinsToPSMs.size(); k++)
    {
      const Set & s = proteinsToPSMs.associations[k];
      proteinsToPSMs.weights[k] = Vector(PSMsToProteins.weights[ s ]).sum();
    }
}

