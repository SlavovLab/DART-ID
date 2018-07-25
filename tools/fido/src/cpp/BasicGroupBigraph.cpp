#include "BasicGroupBigraph.h"

bool gbUseProteinGroupLevelInference = false;

void BasicGroupBigraph::printProteinWeights() const
{
  Array<double> sorted = proteinsToPSMs.weights;
  Array<int> indices = sorted.sort();

  for (int k=0; k<proteinsToPSMs.size(); k++)
    {
      cout << sorted[k] << " " << groupProtNames[ indices[k] ] << endl;
    }
}

void BasicGroupBigraph::trivialGroupProteins()
{
  groupProtNames = Array<Array<string> >(proteinsToPSMs.size());
  originalN = Array<Counter>(proteinsToPSMs.size());

  for (int k=0; k<proteinsToPSMs.size(); k++)
    {
      groupProtNames[k] = Array<string>(1, proteinsToPSMs.names[k]);
      originalN[k] = Counter(1);
    }
}

void BasicGroupBigraph::groupProteinsBy(const Array<Set> & groups)
{
  // remake the list of names and then remake the graph with them collapsed
  groupProtNames = Array<Array<string> > (groups.size());
  originalN = Array<Counter> (groups.size());

  int k;
  for (k=0; k<groups.size(); k++)
    {
      groupProtNames[k] = proteinsToPSMs.names[ groups[k] ];

      // all possible sets of proteins can be present for a group
      if(!gbUseProteinGroupLevelInference){ 
              originalN[k] = Counter( groups[k].size() );
      }else{
      // each group is either present or absent
              originalN[k] = Counter( 1 );
      }
    }

  // remove all but the first of each group from the graph
  for (k=0; k<groups.size(); k++)
    {
      Set reps = groups[k].without( Set::SingletonSet(groups[k][0]) );
      for (Set::Iterator iter = reps.begin(); iter != reps.end(); iter++)
	{
	  disconnectProtein(*iter);
	}
    }

  reindex();
}

void BasicGroupBigraph::groupProteins()
{
  Array<Set> groups = ReplicateIndexer<Set>::replicates(Set::sumSetElements, proteinsToPSMs.associations );

  groupProteinsBy(groups);
}

int BasicGroupBigraph::numberAssociatedProteins(int indexEpsilon) const
{
  int tot = 0;

  const Set & s = PSMsToProteins.associations[indexEpsilon];
  for (int k=0; k<s.size(); k++)
    {
      tot += originalN[ s[k] ].size;
    }
  
  return tot;
}

int BasicGroupBigraph::numberActiveAssociatedProteins(int indexEpsilon, const Array<Counter> & n) const
{
  int tot = 0;

  const Set & s = PSMsToProteins.associations[indexEpsilon];
  for (int k=0; k<s.size(); k++)
    {
      tot += n[ s[k] ].state;
    }
  
  return tot;
}

Array<double> BasicGroupBigraph::probabilityE(const Model & m) const
{
  Array<double> result( PSMsToProteins.size() );

  for (int k=0; k<result.size(); k++)
    {
      result[k] = probabilityEEpsilon(m, k);
    }

  return result;
}

double BasicGroupBigraph::probabilityEEpsilon(const Model & m, int indexEpsilon) const
{
  double tot = 0.0;
  int a = numberAssociatedProteins(indexEpsilon);
  
  for (int k=0; k <= a; k++)
    {
      tot += probabilityEEpsilonGivenActiveAssociatedProteins(m, indexEpsilon, k) * probabilityNumberAssociatedProteins(m, a, k);
    }

  return tot;
}

double BasicGroupBigraph::probabilityEEpsilonGivenActiveAssociatedProteins(const Model & m, int indexEpsilon, int active) const
{
  return 1 - m.probabilityNoEmissionFrom(indexEpsilon, active);
}

double BasicGroupBigraph::probabilityNumberAssociatedProteins(const Model & m, int total, int active) const
{
  return m.probabilityProteins(total, active);
}

double BasicGroupBigraph::probabilityEEpsilonGivenN(const Model & m, int indexEpsilon, const Array<Counter> & n) const
{
  int active = numberActiveAssociatedProteins(indexEpsilon, n);

  return probabilityEEpsilonGivenActiveAssociatedProteins(m, indexEpsilon, active);
}

double BasicGroupBigraph::likelihoodProdTerm(int indexEpsilon, const Model & m, const Array<Counter> & n) const
{
  double probEGivenD = PSMsToProteins.weights[indexEpsilon];
  double probEGivenN = probabilityEEpsilonGivenN(m, indexEpsilon, n);

  double probE = PeptideProphetPriorAtChargeState[ PSMsToProteins.chargeStates[indexEpsilon] ];
  //  cerr << "using prior " << probE << "for c = " << PSMsToProteins.chargeStates[indexEpsilon] << endl;

  double termE = probEGivenD / probE * probEGivenN;
  double termNotE = (1-probEGivenD) / (1-probE) * (1-probEGivenN);
  double term = termE + termNotE;

  return term;
}

double BasicGroupBigraph::likelihoodNGivenD(const Model & m, const Array<Counter> & n) const
{
  double prod = 1.0;

  //  cout << "\tL(N=n | D) for n = " << n << " is ";
  for (int k=0; k<PSMsToProteins.size(); k++)
    {      
      prod *= likelihoodProdTerm(k, m, n);
    }

  //  cout << prod << endl;

  return prod;
}

double BasicGroupBigraph::logLikelihoodNGivenD(const Model & m, const Array<Counter> & n) const
{
  double logProd = 0.0;

  //  cout << "\tL(N=n | D) for n = " << n << " is ";
  for (int k=0; k<PSMsToProteins.size(); k++)
    {      
      logProd += log2( likelihoodProdTerm(k, m, n) );
    }

  //  cout << prod << endl;

  return logProd;
}

double BasicGroupBigraph::probabilityN(const Model & m, const Array<Counter> & n) const
{
  double prod = 1.0;

  for (int k=0; k<n.size(); k++)
    {
      prod *= probabilityNNu(m, n[k]);
    }

  return prod;
}

double BasicGroupBigraph::logNumberOfConfigurations() const
{
  double result = 0.0;

  //  cout << "\tGetting logNumber configs for " << originalN << endl;
  // num = prod(originalN[k]*#peptides,k)
  for (int k=0; k<originalN.size(); k++)
    {
      result += log2(originalN[k].size+1);
    }

  //  cout << "\t\twas = " << result << endl;

  //  return result + log2(PSMsToProteins.size());
  return result;
}

double BasicGroupBigraph::probabilityNNu(const Model & m, const Counter & nNu) const
{
  return m.probabilityProteins(nNu.size, nNu.state);
}

double BasicGroupBigraph::probabilityNGivenD(const Model & m, const Array<Counter> & n) const
{
  //  return likelihoodNGivenD(m, n) * probabilityN(m, n) / likelihoodConstant(m);

  // using cached functor
  //  cout << "Calling functor" << endl;
  //  cout << likelihoodConstantCachedFunctor.name << endl;
  //  cout << "correct object address is " << this << endl;

  // working version, but numerically unstable
  //  double like = likelihoodNGivenD(m, n) * probabilityN(m, n) / likelihoodConstantCachedFunctor(m, this);

  // log version
  double logLike= logLikelihoodNGivenD(m,n) + log2(probabilityN(m,n)) - logLikelihoodConstantCachedFunctor(m,this);

  /***
  cout << "L = " << like << endl;
  cout << "2^logL = " << pow(2.0, logLike) << endl << endl;

  cout << "L(N|D) = " << likelihoodNGivenD(m,n) << endl;
  cout << "2^logL(N|D) = " << pow( 2.0, logLikelihoodNGivenD(m,n) ) << endl << endl;
  ***/

  //return like;
  return pow(2.0, logLike);
}

double BasicGroupBigraph::logLikelihoodConstant(const Model & m) const
{
  double result = 0.0;
  bool starting = true;

  Array<Counter> n = originalN;
  for (Counter::start(n); Counter::inRange(n); Counter::advance(n))
    {
      double L = logLikelihoodNGivenD(m, n);
      double p = log2(probabilityN(m, n));
      double logLikeTerm = L+p;

      if ( starting )
	{
	  starting = false;
	  result = logLikeTerm;
	}
      else
	{
	  result = Numerical::logAdd(result, logLikeTerm);
	}
    }

  return result;
}

double BasicGroupBigraph::likelihoodConstant(const Model & m) const
{
  double result = 0.0;

  Array<Counter> n = originalN;
  for (Counter::start(n); Counter::inRange(n); Counter::advance(n))
    {
      double L = likelihoodNGivenD(m, n);
      double p = probabilityN(m, n);
      double term = L*p;
      //      cout << "\tL and p are: " << L << ", " << p << endl;
      //      cout << "\t\tterm is " << term << endl;

      //      cout << "\tif you'd used logL: " << logLikelihoodNGivenD(m,n) << endl;
      result += term;
    }

  //  if ( isinf( result ) )
  //    print();
  //    displayDotty("Checker");

  //  cout << "\tWith result " << result << endl;
  return result;
}

Array<double> BasicGroupBigraph::probabilityXGivenD(bool xState, const Model & m)
{
  Array<Counter> n = originalN;

  Vector result;

  for (Counter::start(n); Counter::inRange(n); Counter::advance(n))
    {
      Vector term = probabilityNGivenD(m, n) * Vector( xCorrection(xState, n) );

      if ( result.size() == 0 )
	{
	  result = term;
	}
      else
	{
	  result += term;
	}
    }

  return result.unpack();
}

Array<double> BasicGroupBigraph::probabilityXAndEGivenD(bool xState, const Model & m)
{
  Array<Counter> n = originalN;

  Vector result;

  for (Counter::start(n); Counter::inRange(n); Counter::advance(n))
    {
      Vector term = probabilityNGivenD(m, n) * Vector( xAndECorrection(xState, m, n) );

      if ( result.size() == 0 )
	{
	  result = term;
	}
      else
	{
	  result += term;
	}
    }

  return result.unpack();
}

Array<double> BasicGroupBigraph::xCorrection(bool xState, const Array<Counter> & n)
{
  Array<double> result(PSMsToProteins.size());
  for (int k=0; k<result.size(); k++)
    {
      result[k] = xCorrectionEpsilon(k, xState, n);
    }

  return result;
}

double BasicGroupBigraph::xCorrectionEpsilon(int indexEpsilon, bool xState, const Array<Counter> & n)
{
  const Set & s = PSMsToProteins.associations[indexEpsilon];
  for (Set::Iterator iter = s.begin(); iter != s.end(); iter++)
    {
      if ( n[*iter].state > 0 )
	{
	  // when xState is true, return 1
	  // when xState is false, return 0
	  return double(xState);
	}
    }

  // when xState is true, return 0
  // when xState is false, return 1
  return 1.0 - double(xState);
}

Array<double> BasicGroupBigraph::probabilityEGivenD(const Model & m)
{
  Array<Counter> n = originalN;

  //  cout << "originalN = " << originalN << endl;

  Vector result;

  for (Counter::start(n); Counter::inRange(n); Counter::advance(n))
    {
      Vector term = probabilityNGivenD(m, n) * Vector( eCorrection(m, n) );

      if ( result.size() == 0 )
	{
	  result = term;
	}
      else
	{
	  result += term;
	}
    }

  return result.unpack();
}

Array<double> BasicGroupBigraph::estimateAlphaGivenD(const Model & m)
{
  getProteinProbs(m);

  Array<Counter> n = originalN;
  Vector result;

  for (Counter::start(n); Counter::inRange(n); Counter::advance(n))
    {
      //      cout << "adding term with " << (probabilityNGivenD(m, n) * Vector(alphaEstimateCorrection(m, n)) ).unpack() << endl;

      // note: the alphaEstimateCorrection values are numerically
      // large, but the product between them and the coefficient
      // probabilityNGivenD is small; therefore, it may be a good idea
      // to log both, add them, and then exponentiate to create each
      // term
      Vector term = probabilityNGivenD(m, n) * Vector( alphaEstimateCorrection(m, n) );

      if ( result.size() == 0 )
	{
	  result = term;
	}
      else
	{
	  result += term;
	}
    }

  return result.unpack();
}

Array<double> BasicGroupBigraph::alphaEstimateCorrection(const Model & m, const Array<Counter> & n) const
{
  Array<double> result(PSMsToProteins.size());

  for (int k=0; k<result.size(); k++)
    {
      result[k] = alphaEstimateCorrectionEpsilon(k, m, n);
    }

  return result;
}

double BasicGroupBigraph::alphaEstimateCorrectionEpsilon(int indexEpsilon, const Model & m, const Array<Counter> & n) const
{
  double sum = 0.0;

  const Set & s = PSMsToProteins.associations[indexEpsilon];
  for (Set::Iterator iter = s.begin(); iter != s.end(); iter++)
    {
      double term = alphaEstimateCorrectionNuEpsilon( *iter, indexEpsilon, m, n, probabilityR[*iter] );
      sum += term * n[*iter].size;
    }

  return sum / (numberAssociatedProteins(indexEpsilon));
}

double BasicGroupBigraph::alphaEstimateCorrectionNuEpsilon(int nu, int indexEpsilon, const Model & m, const Array<Counter> & n, double probRRhoGivenD) const
{
  double probEGivenD = PSMsToProteins.weights[indexEpsilon];
  double probE = PeptideProphetPriorAtChargeState[ PSMsToProteins.chargeStates[indexEpsilon] ];
  //  cerr << "using prior " << probE << "for c = " << PSMsToProteins.chargeStates[indexEpsilon] << endl;

  double termE = probEGivenD / probE;

  double likeCorrection = termE / likelihoodProdTerm(indexEpsilon, m, n);
  double rCorrection = probabilityRRhoGivenN(nu, n);

  return likeCorrection * rCorrection * m.associatedEmission(indexEpsilon) / probRRhoGivenD;
}

Array<double> BasicGroupBigraph::xAndECorrection(bool xState, const Model & m, const Array<Counter> & n)
{
  Array<double> result(PSMsToProteins.size());

  for (int k=0; k<result.size(); k++)
    {
      result[k] = xCorrectionEpsilon(k, xState, n) * eCorrectionEpsilon(k, m, n);
    }

  return result;
}

Array<double> BasicGroupBigraph::eCorrection(const Model & m, const Array<Counter> & n)
{
  Array<double> result(PSMsToProteins.size());

  for (int k=0; k<result.size(); k++)
    {
      result[k] = eCorrectionEpsilon(k, m, n);
    }

  return result;
}

double BasicGroupBigraph::eCorrectionEpsilon(int indexEpsilon, const Model & m, const Array<Counter> & n)
{
  double probEGivenD = PSMsToProteins.weights[indexEpsilon];
  double probEGivenN = probabilityEEpsilonGivenN(m, indexEpsilon, n);

  double probE = PeptideProphetPriorAtChargeState[ PSMsToProteins.chargeStates[indexEpsilon] ];
  //  cerr << "using prior " << probE << "for c = " << PSMsToProteins.chargeStates[indexEpsilon] << endl;

  double termE = probEGivenD / probE * probEGivenN;
  double termNotE = (1-probEGivenD) / (1-probE) * (1-probEGivenN);
  double term = termE + termNotE;

  return termE / term;
}

Array<double> BasicGroupBigraph::MAPRGivenD(const Model & m)
{
  Array<Counter> n = originalN;

  Array<Counter> best;
  double bestValue = -1;

  for (Counter::start(n); Counter::inRange(n); Counter::advance(n))
    {
      double term = likelihoodNGivenD(m, n) * probabilityN(m, n);

      //      cout << "Current score: " << term << endl;
      //      cout << proteinGroupNames() << endl;
      //      cout << '\t' << n << endl;

      if ( term > bestValue)
	{
	  bestValue = term;
	  best = n;
	}
    }

  //  cout << "best was " << bestValue << " from " << best << endl;

  Array<double> result(n.size());
  
  for (int k=0; k<n.size(); k++)
    {
      result[k] = double(best[k].state > 0);
    }
  
  return result;
}

Array<double> BasicGroupBigraph::probabilityRGivenD(const Model & m)
{
  Array<Counter> n = originalN;

  //  cout << "originalN = " << originalN << endl;

  Vector result;

  for (Counter::start(n); Counter::inRange(n); Counter::advance(n))
    {
      /***
      cout << "Current iter: " << endl;
      cout << '\t' << n << endl;
      cout << '\t' << probabilityNGivenD(m,n);
      cout << '\t' << probabilityRGivenN(n) << endl << endl;
      ***/

      Vector term = probabilityNGivenD(m, n) * Vector( probabilityRGivenN(n) );

      if ( result.size() == 0 )
	{
	  result = term;
	}
      else
	{
	  result += term;
	}
    }

  return result.unpack();
}

Array<double> BasicGroupBigraph::probabilityRGivenN(const Array<Counter> & n)
{
  Array<double> result(n.size());

  for (int k=0; k<result.size(); k++)
    {
      result[k] = probabilityRRhoGivenN(k, n);
    }

  return result;
}

double BasicGroupBigraph::probabilityRRhoGivenN(int indexRho, const Array<Counter> & n) const
{
  const Counter & c = n[indexRho];

  return double(c.state) / c.size;
}

void BasicGroupBigraph::getProteinProbs(const Model & m)
{
  //  cout << m << endl;
  probabilityR = probabilityRGivenD(m);
  //  cout << "BasicGroupBigraph::probR = " << probabilityR << endl;
}

void BasicGroupBigraph::getMAPProteins(const Model & m)
{
  //  cout << m << endl;
  probabilityR = MAPRGivenD(m);
  //  cout << "BasicGroupBigraph::probR = " << probabilityR << endl;
}

double BasicGroupBigraph::probabilityEEpsilonOverAllAlphaBeta(const GridModel & gm, int indexEpsilon) const
{
  GridModel localModel( gm );

  double val = 0.0;
  for ( localModel.start(); localModel.inRange(); localModel.advance() )
    {
      val += probabilityEEpsilon(localModel, indexEpsilon);
    }

  val /= localModel.getCount();

  return val;
}

Array<double> BasicGroupBigraph::probabilityEOverAllAlphaBeta(const GridModel & gm) const
{
  Array<double> result( PSMsToProteins.size() );

  for ( int k=0; k<result.size(); k++)
    {
      result[k] = probabilityEEpsilonOverAllAlphaBeta(gm, k);
    }

  return result;
}

/*
double BasicGroupBigraph::likelihoodAlphaBetaGivenD(const GridModel & gm) const
{
  return likelihoodConstantCachedFunctor(gm, this);
}
*/

double BasicGroupBigraph::logLikelihoodAlphaBetaGivenD(const GridModel & gm) const
{
  return logLikelihoodConstantCachedFunctor(gm, this);
}

