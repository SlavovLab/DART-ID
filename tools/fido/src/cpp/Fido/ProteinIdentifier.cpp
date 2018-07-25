#include "ProteinIdentifier.h"

Array<double> ProteinIdentifier::PeptideProphetPriorAtChargeState = Array<double>(10, 0.0);

ProteinIdentifier::ProteinIdentifier()
{
  //  ProteinThreshold = 1e-2;
  ProteinThreshold = 1e-2;
  
  //  PeptideThreshold = 1e-9;
  //  PeptideThreshold = 9e-3;
  //  PeptideThreshold = 1e-9;

  PeptideThreshold = 1e-3;

  PeptideProphetPriorAtChargeState[0] = 0.10; // when there is no information
}

