#include "PeptideFeatureReader.h"

PeptideFeatureReader::PeptideFeatureReader(double g):
  //  sparsify(1e-2)
  sparsify( -log(0.75)/log(2.718281728459)*g )
{
  rbfGamma = g;
}

void PeptideFeatureReader::read(char*fname)
{
  ifstream fin(fname);
  string buff, wd;

  getline(fin, buff);
  istringstream featureIst(buff);

  cout << "feature line was: " << buff << endl;

  // dump the first feature name it should be the peptide feature
  // (sequence of the peptide)
  featureIst >> wd;
  while (featureIst >> wd)
    {
      featureNames.add(wd);
    }

  double value;
  while ( getline(fin, buff) )
    {
      istringstream ist(buff);
      ist >> wd;

      Array<double> featureArray;
      while (ist >> value)
	{
	  featureArray.add(value);
	}
       
      if ( featureArray.size() != featureNames.size() )
	{
	  cerr << "PeptideFeatureReader error: length of line doesn't match length of feature names" << endl;
	  cerr << "\tline was: " << buff << endl << endl;
	  exit(1);
	}

      if ( sequences.size() > 0 && sequences.back() >= wd )
	{
	  cerr << "PeptideFeatureReader error: sequences read by PeptideFeatureReader must be sorted (ascending) and unique" << endl << endl;
	  exit(1);
	}

      sequences.add(wd);
      featureVectors.add( Vector(featureArray) );
   }
}

Array<int> PeptideFeatureReader::findPeptideSequences(const Array<string> & pepSeqs)
{
  Array<int> result(pepSeqs.size());

  for (int k=0; k<pepSeqs.size(); k++)
    {
      int loc = sequences.sortedFind(pepSeqs[k]);

      if (loc == -1)
	{
	  cerr << "Error: peptide " << pepSeqs[k] << " was not found in the pepFtrs file" << endl;
	  exit(1);
	}

      result[k] = loc;
    }

  return result;
}

void PeptideFeatureReader::storeProbabilities(const Array<string> & pepSeqs, const Array<double> & xAndEGivenD, const Array<double> & xGivenD)
{
  Array<int> locs = findPeptideSequences(pepSeqs);
  probXAndEGivenD = Array<Array<double> >(sequences.size());
  probXGivenD = Array<Array<double> >(sequences.size());

  usedIndices.clear();

  for (int k=0; k<locs.size(); k++)
    {
      int i = locs[k];
      probXAndEGivenD[ i ].add(xAndEGivenD[k]);
      probXGivenD[ i ].add(xGivenD[k]);

      usedIndices.insert(i);
    }
}

Array<double> PeptideFeatureReader::estimateAlphasSlow(double meanAlpha, const Array<string> & pepSeqs, const Array<double> & xAndEGivenD, const Array<double> & xGivenD)
{
  storeProbabilities(pepSeqs, xAndEGivenD, xGivenD);

  Array<double> all_answers(sequences.size(), -1.0);

  Array<int> locs = findPeptideSequences(pepSeqs);

  // cache the total numerator and denominator scores for each group of identical peptides
  cerr << "caching..." << endl;
  Array<double> cachedMultiNumers(sequences.size(), -1.0), cachedMultiDenoms(sequences.size(), -1.0);
  for (set<int>::const_iterator iter=usedIndices.begin(); iter != usedIndices.end(); iter++)
    {
      int u = *iter;

      double numer = 0.0;
      double denom = 0.0;

      //      cerr << "Caching " << "index " << u << " " << sequences[u] << endl;

      for (int j=0; j<probXAndEGivenD[u].size(); j++)
	{
	  numer += probXAndEGivenD[ u ][j];
	  denom += probXGivenD[ u ][j];

	  //	  cerr << "\t" << "identical to " << numer << "/" << denom << endl;
	}

      cachedMultiNumers[u] = numer;
      cachedMultiDenoms[u] = denom;

      //      cerr << "\tresults in " << numer << "/" << denom << endl;

      //      cout << "\t\t" << numer << " / " << denom << endl;
    }

  cerr << "smoothing the alphas..." << endl;
  for (set<int>::const_iterator iter_k=usedIndices.begin(); iter_k != usedIndices.end(); iter_k++)
    {
      int k = *iter_k;
      double numer = 0.0;
      double denom = 0.0;

      //      cerr << "Smoothing " << sequences[k] << endl;
      for (set<int>::const_iterator iter_i=usedIndices.begin(); iter_i != usedIndices.end(); iter_i++)
	{
	  int i = *iter_i;
	  double Kki = K(k, i);

	  numer += cachedMultiNumers[i] * Kki;
	  denom += cachedMultiDenoms[i] * Kki;

	  if ( Kki > 0.05 && denom > 0.05 )
	    {
	      //	      cerr << "\t" << Kki << " x " << cachedMultiNumers[i] << "/" << cachedMultiDenoms[i] << " from index " << i << " " << sequences[i] << endl;
	    }
	}
      //      cerr << "\tresults in " << numer << "/" << denom << endl;

      all_answers[k] = computeProb(meanAlpha, numer, denom);
      //      cout << "\t\t" << all_answers[k] << endl;
      //      cout << "\t\t\t" << numer << " / " << denom << ", " << meanAlpha << endl;
      //      cout << sequences[k] << " = " << numer << "/" << denom << " = " << result[k] << endl;
    }
  
  cerr << "indexing them for results..." << endl;
  Array<double> result(pepSeqs.size());
  for (int k=0; k<result.size(); k++)
    {
      result[k] = all_answers[ locs[k] ];
    }

  return result;
}

void PeptideFeatureReader::slowComputePeptideSimilarity()
{
  // note: this is currently a quadratic (in number of peptides)
  // memory requirement

  Array<Array<double> > similarityArrayMat(featureVectors.size());
  for (int i=0; i<featureVectors.size(); i++)
    {
      similarityArrayMat[i] = Array<double>(featureVectors.size());
      int j;
      for (j=0; j<=i; j++)
	{
	  similarityArrayMat[i][j] = K(i,j);
	  similarityArrayMat[j][i] = similarityArrayMat[i][j];
	}
    }

  similarity = Matrix(similarityArrayMat);
}

Array<double> PeptideFeatureReader::estimateAlphas(double meanAlpha, const Array<double> & probabilityXAndEGivenD, const Array<double> & probabilityXGivenD)
{
  /***
  Set nonzero;

  int i;
  cout << "checking sizes are ==: " << probabilityXAndEGivenD.size() << " " << probabilityXGivenD.size() << " " << featureVectors.size() << endl;

  for (i=0; i<probabilityXGivenD.size(); i++)
    {
      if ( sparsify.isNonzero(probabilityXGivenD[i]) )
	nonzero.add(i);
    }

  cout << "Went from " << i-1 << " cols --> " << nonzero.size() << " cols" << endl;
  cout << "Performing actual alpha computation" << endl;
  ***/

  Array<double> result(probabilityXAndEGivenD.size());

  for (int k=0; k<result.size(); k++)
    {
      const Array<int> & s = similarPeptides[k];
      //      cout << "links for " << k << " are " << s << endl;
      //      const Array<int> & s = ;

      double numer = 0.0;
      double denom = 0.0;
      for (int i=0; i<s.size(); i++)
	{

	  double Ki = K(k, s[i]);

	  //	  cout << "\tnumer adds " << probabilityXAndEGivenD[ s[i] ] * Ki << endl;
	  //	  cout << "\tdenom adds " << probabilityXGivenD[ s[i] ] * Ki << endl;

	  numer += probabilityXAndEGivenD[ s[i] ] * Ki;
	  denom += probabilityXGivenD[ s[i] ] * Ki;
	}

      result[k] = computeProb(meanAlpha, numer, denom);
      //      cout << sequences[k] << " = " << numer << "/" << denom << " = " << result[k] << endl;
    }

  return result;
}

double PeptideFeatureReader::computeProb(double meanAlpha, double positives, double trials) const
{
  double est = positives/trials;
  return est;

  // note: the following code is not numerically stable when using a
  // large number of positives and trials (ie when you're smoothing
  // aggressively)
  double propLikeEst = pow(est, positives)*pow(1-est, trials-positives);
  double propLikeMean = pow(meanAlpha, positives)*pow(1-meanAlpha, trials-positives);

  double prior = priorProb(trials);

  double alphaStar = (propLikeEst*est*prior + propLikeMean*meanAlpha*(1-prior))/(propLikeEst*prior+propLikeMean*(1-prior));

  return alphaStar;
  //  return min(max( est, .01 ), .6 );
}

double PeptideFeatureReader::priorProb(double trials) const
{
  // slowness is a constant in (0,1]; if it is close to 0, the prior
  // turns to 1.0 almost instantly. If it is close to 1, the prior
  // changes very slowly

  double slowness = .8;
  return 1-pow(slowness, trials);
}

double PeptideFeatureReader::K(int i, int j) const
{
  //  return exp( pow((featureVectors[i] - featureVectors[j]).norm(), 2.0) / rbfGamma );
  return exp( -diffNormSquared(featureVectors[i], featureVectors[j]) / rbfGamma );
}

void PeptideFeatureReader::normalizePercentiles()
{
  Array<Array<double> > featureArrays(featureVectors.size(), Array<double>(featureNames.size()) );

  int k,j;
  for (k=0; k<featureNames.size(); k++)
    {
      Array<double> colK = percentileFeature(k);

      for (j=0; j<featureVectors.size(); j++)
	{
	  featureArrays[j][k] = colK[j];
	}
    }

  for (j=0; j<featureVectors.size(); j++)
    {
      featureVectors[j] = Vector(featureArrays[j]);
    }
}

Array<double> PeptideFeatureReader::percentileFeature(int feature)
{
  Array<double> fA( featureVectors.size() );
  
  for (int k=0; k<fA.size(); k++)
    {
      fA[k] = featureVectors[k][feature];
    }

  Array<double> sorted = fA;
  Array<int> order = sorted.sortA();

  Array<double> sortedPercentiles = getNormalizedRanking(sorted);
  
  // from now on, fA is only used as the result
  for (int k=0; k<sorted.size(); k++)
    {
      fA[ order[k] ] = sortedPercentiles[ k ];
    }

  return fA;
}

/***
void PeptideFeatureReader::linkSimilar()
{
  KDMethod splitter(sparsify, orderedFeatureVectorRefs);
  similarPeptides = splitter.getSimilar();
}
***/
