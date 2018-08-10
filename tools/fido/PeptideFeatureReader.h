#ifndef _PEPTIDEFEATUREREADER_H
#define _PEPTIDEFEATUREREADER_H

#include "Matrix.h"
#include <fstream>
#include <math.h>
#include <sstream>
#include "Random.h"
#include "KDMethod.h"
#include <set>

class PeptideFeatureReader
{
 public:
  PeptideFeatureReader(double g);
  void read(char*fname);
  void reorder(const Array<string> & pepNames);
  void slowComputePeptideSimilarity();
  double K(int i, int j) const;

  void normalizePercentiles();
  Array<double> percentileFeature(int feature);

  void linkSimilar();

  Array<int> findPeptideSequences(const Array<string> & pepSeqs);
  void storeProbabilities(const Array<string> & pepSeqs, const Array<double> & xAndEGivenD, const Array<double> & xGivenD);
  Array<double> estimateAlphasSlow(double meanAlpha, const Array<string> & pepSeqs, const Array<double> & xAndEGivenD, const Array<double> & xGivenD);

  Array<double> estimateAlphas(double meanAlpha, const Array<double> & probabilityXAndEGivenD, const Array<double> & probabilityXGivenD);
  double computeProb(double meanAlpha, double positives, double trials) const;
  double priorProb(double trials) const;

  int findRankBound(const Array<double> & sorted, int k1)
  {
    int k2;
    for (k2 = k1; k2<sorted.size()-1; k2++)
      {
	if ( sorted[k2+1] != sorted[k1] )
	  break;
      }

    // there are two cases arriving here:
    // 1. break statement (trivially OK)
    // 2. k2 = a.size()-1 and all between k1, k2 are equal

    return k2;
  }

  Array<double> getNormalizedRanking(const Array<double> & sorted)
    {
      Array<double> result(sorted.size(), -1);

      double window = 0;
      for (int k1 = 0; k1 < sorted.size(); )
	{
	  int k2 = findRankBound(sorted, k1);

	  double windowSize = k2-k1 + 1;
	  double averageRank = ( ( window + windowSize )*( window + windowSize + 1 ) - window * (window + 1) ) / ( 2*windowSize );
	  //      cout << "Giving window " << k1 << ", " << k2 << "\t" << averageRank << endl;

	  for (int k=k1; k<=k2; k++)
	    {
	      result[k] = averageRank / sorted.size();
	    }
      
	  window += k2 - k1 + 1;
	  k1 = k2+1;
	}

      return result;
    }

  Array<string> featureNames;
  Array<string> sequences;
  Array<Vector> featureVectors;
  Array<Array<double> > probXAndEGivenD, probXGivenD;
  set<int> usedIndices;

  Array<Vector> percentileFeatures;
  Array<Array<int> > similarPeptides;

  Matrix similarity;
  double rbfGamma;
  Numerical sparsify;
};


#endif
