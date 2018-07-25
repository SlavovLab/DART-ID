#include "EMBigraph.h"

void EMBigraph::EM(int numberOfIterations)
{
  //  cout << "Starting EM" << endl;
  //  cout << "Found " << proteinsToPSMs.size() << " proteins" << endl;

  int k;
  for (k=0; k<proteinsToPSMs.size(); k++)
    proteinsToPSMs.weights[k] = 1.0;

  for (k = 0; k<numberOfIterations; k++)
    EMIteration();
}

void EMBigraph::EMIteration()
{
  //  cout << "\tEM Iter" << endl;
  Array<double> newProteinWeights( proteinsToPSMs.size(), 0.0 );

  int k;
  for (k=0; k<newProteinWeights.size(); k++)
    {
      newProteinWeights[k] = probR(k);
      //      proteinsToPSMs.weights[k] = probR(k);
    }

  proteinsToPSMs.weights = newProteinWeights;
}

double EMBigraph::probR(int index)
{
  double prod = 1.0;

  const Set & s = proteinsToPSMs.associations[index];
  for (int k=0; k<s.size(); k++)
    {
      int eIndex = s[k];
      double eWeight = PSMsToProteins.weights[eIndex];
      double term = (1-eWeight) + eWeight * existsOtherRForE(index, eIndex);

      //      prod *= term;
      prod *= pow(term, expectedCountE(eIndex) );
    }

  return 1 - prod;
}

double EMBigraph::expectedCountE(int index)
{
  return Vector(proteinsToPSMs.weights[ PSMsToProteins.associations[index] ] ).sum();
}

double EMBigraph::existsOtherRForE(int rIndex, int eIndex)
{
  double prod = 1.0;

  const Set & s = PSMsToProteins.associations[eIndex];
  for (int k=0; k<s.size(); k++)
    {
      int protIndex = s[k];

      if ( protIndex != rIndex )
	prod *= proteinsToPSMs.weights[protIndex];
    } 

  return 1 - prod;
}


