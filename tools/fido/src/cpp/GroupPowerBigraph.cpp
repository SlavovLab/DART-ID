#include "GroupPowerBigraph.h"

double GroupPowerBigraph::LOG_MAX_ALLOWED_CONFIGURATIONS = 18;
//double GroupPowerBigraph::LOG_MAX_ALLOWED_CONFIGURATIONS = 1500000000;

Array<double> GroupPowerBigraph::proteinProbs()
{
  Array<double> result;
  
  for (int k=0; k<subgraphs.size(); k++)
    {
      //cout << "Working on subgraph " << k << endl;
      //      cout << "using model " << subgraphModels[k].alphaArray << " " << subgraphModels[k] << endl;

      //      if ( subgraphs[k].proteinsToPSMs.size() > 1 )
      //      cerr << "Working on subgraph with " << subgraphs[k].proteinsToPSMs.size() << " proteins (" << subgraphs[k].proteinsToPSMs.names << ") and " << subgraphs[k].PSMsToProteins.size() << " peptides"<< endl;
      //      cerr << subgraphs[k].PSMsToProteins.weights << endl << endl;

      subgraphs[k].getProteinProbs(subgraphModels[k]);
      result.append( subgraphs[k].proteinProbabilities() );
    }

  // clip into region [0,1] to correct numerical errors
  for (int j=0; j<result.size(); ++j)
    result[j] = max(0.0, min(1.0, result[j]));

  return result;
}

Array<double> GroupPowerBigraph::MAPProteins()
{
  Array<double> result;
  
  for (int k=0; k<subgraphs.size(); k++)
    {
      //      cout << "Working on subgraph " << k << endl;
      //      cout << "using model " << subgraphModels[k].alphaArray << " " << subgraphModels[k] << endl;

      //      if ( subgraphs[k].proteinsToPSMs.size() > 5 )
      //	cerr << "Working on subgraph with " << subgraphs[k].proteinsToPSMs.size() << " proteins and " << subgraphs[k].PSMsToProteins.size() << " peptides"<< endl;

      subgraphs[k].getMAPProteins(subgraphModels[k]);
      result.append( subgraphs[k].proteinProbabilities() );
    }

  return result;
}

void GroupPowerBigraph::printGraphs()
{
  for (int k=0; k<subgraphs.size(); k++)
    {
      cout << "Subgraph " << k << endl;

      subgraphs[k].printGraph();

      cout << endl << endl;
    }
}

void GroupPowerBigraph::printPeptideDetectabilities()
{
  Array<Array<string> > prots;
  Array<string> peps;
  Array<double> pepProph;

  int k;
  for (k=0; k<subgraphs.size(); k++)
    {
      for (int j=0; j<subgraphs[k].PSMsToProteins.size(); j++)
	{
	  prots.add( subgraphs[k].proteinsToPSMs.names[ subgraphs[k].PSMsToProteins.associations[j] ] );
	  pepProph.add( subgraphs[k].PSMsToProteins.weights[j] );
	}
      peps.append( subgraphs[k].PSMsToProteins.names );
    }

  // should maybe use prEGivenD, but that won't force emission rather
  // than creation from noise
  Array<double> xAndE = probabilityXAndEGivenD();
  Array<double> x = probabilityXGivenD();

  for (k=0; k<peps.size(); k++)
    {
      //      cerr << peps[k] << '\t' << prots[k] << '\t' << xAndE[k] << '\t' << x[k] << endl;
      cout << peps[k] << '\t' << xAndE[k] << '\t' << x[k] << endl;
    }
}

Array<double> GroupPowerBigraph::probabilityXAndEGivenD()
{
  Array<double> result;

  for (int k=0; k<subgraphs.size(); k++)
    {
      result.append( subgraphs[k].probabilityXAndEGivenD( true , subgraphModels[k] ) );
    }

  return result;
}

Array<double> GroupPowerBigraph::probabilityXGivenD()
{
  Array<double> result;

  for (int k=0; k<subgraphs.size(); k++)
    {
      result.append( subgraphs[k].probabilityXGivenD( true , subgraphModels[k] ) );
    }

  return result;
}

void GroupPowerBigraph::refreshModels()
{
  cout << "estimating alphas" << endl;
  Array<double> probXAndEGivenD = probabilityXAndEGivenD();
  Array<double> probXGivenD = probabilityXGivenD();

  double mean = meanAlpha(probXAndEGivenD, probXGivenD);

  Array<double> alphas = pepFtrRdr.estimateAlphasSlow(mean, peptideNames(), probXAndEGivenD, probXGivenD);

  cout << "Performing other parts of refresh" << endl;
  int i = 0;
  for (int k=0; k<subgraphs.size(); k++)
    {
      Array<double> subgraphAlphas(subgraphs[k].PSMsToProteins.size());
      for (int j=0; j<subgraphs[k].PSMsToProteins.size(); j++, i++)
	{
	  subgraphAlphas[j] = alphas[i];
	}

      if ( subgraphModels[k].alphaArray.size() != subgraphAlphas.size() )
	cerr << "error: alpha sizes don't match" << endl;

      DetectabilityModel & dm = subgraphModels[k];
      dm = DetectabilityModel(subgraphAlphas, dm.beta, dm.gamma);
      
      //      cout << k << ":" << dm.alphaArray << " " << dm << endl;
      subgraphs[k].refreshCache();
    }

  //  showDetectabilities();

  //  cout << "finished refreshing models" << endl;
}

double GroupPowerBigraph::meanAlpha(const Array<double> & probXAndEGivenD, const Array<double> & probXGivenD)
{
  double numer = 0.0, denom = 0.0;

  for (int k=0; k<probXGivenD.size(); k++)
    {
      numer += probXAndEGivenD[k];
      denom += probXGivenD[k];
    }
  
  //  cout << "estimated alpha (overall) is " << numer/denom << endl;
  return numer/denom;
}

void GroupPowerBigraph::getProteinProbs()
{
  probabilityR = proteinProbs();
}

void GroupPowerBigraph::getMAPProteins()
{
  probabilityR = MAPProteins();
}

/***
double GroupPowerBigraph::likelihoodAlphaBetaGivenD(const DetectabilityModel & myGM) const
{
  return pow(2.0, logLikelihoodAlphaBetaGivenD(myGM) );

  double prod = 1.0;

  for (int k=0; k<subgraphs.size(); k++)
    {
      prod *= subgraphs[k].likelihoodAlphaBetaGivenD(myGM);
    }

  return prod / pow( 1-myGM.spontaneousEmission() , numberClones);
}
***/

double GroupPowerBigraph::logLikelihoodAlphaBetaGivenD() const
{
  //  cout << "\tGetting logLike(alpha, beta | D) with " << (Model)myGM << endl;

  double sum = 0.0;

  for (int k=0; k<subgraphs.size(); k++)
    {
      sum += subgraphs[k].logLikelihoodConstantCachedFunctor(subgraphModels[k], & subgraphs[k]);
    }

  return sum - numberClones * log2( 1-subgraphModels[0].spontaneousEmission() );
}

/***
double GroupPowerBigraph::probabilityAlphaBetaGivenD(const DetectabilityModel & myGM) const
{
  // using cached functor
  return pow(2.0, logLikelihoodAlphaBetaGivenD(myGM) - sumLogLikelihoodOverAllAlphaBetaCachedFunctor(myGM, this) );
}

double GroupPowerBigraph::sumLogLikelihoodOverAllAlphaBeta(const DetectabilityModel & myGM) const
{
  DetectabilityModel local = myGM;

  double result = 0.0;
  bool starting = true;

  for (local.start(); local.inRange(); local.advance())
    {
      double logLike = logLikelihoodAlphaBetaGivenD(local);
     
      if ( starting ) 
	{
	  starting = false;
	  result = logLike;
	}
      else
	{
	  result = Numerical::logAdd(result, logLike);
	}
    }

  return result;
}

void GroupPowerBigraph::getProteinProbsOverAllAlphaBeta()
{
  probabilityR = proteinProbsOverAllAlphaBeta();
}


Array<double> GroupPowerBigraph::proteinProbsOverAllAlphaBeta()
{
  Vector cumulative;

  DetectabilityModel local = gm;
  for ( local.start(); local.inRange(); local.advance() )
    {
      double prob = probabilityAlphaBetaGivenD(local);

      // hack for efficiency
      if ( prob > 1e-5 )
	{
	  Array<double> protProbs = proteinProbs( local );
	  Vector posteriorsForCurrentAlphaBeta = prob * Vector( protProbs );

	  if ( cumulative.size() == 0 )
	    cumulative = posteriorsForCurrentAlphaBeta;
	  else
	    cumulative += posteriorsForCurrentAlphaBeta;
	}

    }

  return cumulative.unpack();
}
***/

void GroupPowerBigraph::getGroupProtNames()
{
  int k,j;
  for (k=0; k<subgraphs.size(); k++)
    {
      const BasicGroupBigraph & bgb = subgraphs[k];

      for (j=0; j<bgb.proteinsToPSMs.size(); j++)
	{
	  groupProtNames.add( bgb.proteinGroupNames()[j] );
	}
    }
}

pair<Array<Array<string> >, Array<double> > GroupPowerBigraph::getDescendingProteinsAndWeights() const
{
  Array<double> sorted = probabilityR;
  int k;
  for (k=0; k<sorted.size(); k++)
    {
      if ( isnan( sorted[k] ) )
	{
	  cerr << "error: found nan in GroupPowerBigraph::getDescendingProteinWeights" << endl;
	  cerr << "\tprotein " << groupProtNames[k] << endl;
	}
    }

  Array<int> indices = sorted.sort();

  Array<Array<string> > groupNames;
  Array<double> probabilities;
  for (k=0; k<sorted.size(); k++)
    {
      probabilities.add(sorted[k]);
      groupNames.add(groupProtNames[ indices[k] ] );
    }
  return pair<Array<Array<string> >, Array<double> >(groupNames, probabilities);
}

void GroupPowerBigraph::printProteinWeights() const
{
  //  cerr << "Printing prot weights" << endl;
  Array<double> sorted = probabilityR;
  int k;
  for (k=0; k<sorted.size(); k++)
    {
      if ( isnan( sorted[k] ) )
	{
	  cerr << "error: found nan in GroupPowerBigraph::printProteinWeights" << endl;
	  cerr << "\tprotein " << groupProtNames[k] << endl;
	}
    }

  Array<int> indices = sorted.sort();

  for (k=0; k<sorted.size(); k++)
    {
      cout << sorted[k] << " " << groupProtNames[ indices[k] ] << endl;
    }

  cout << "0.0 " << severedProteins << endl;
}

void GroupPowerBigraph::readFromMCMC(istream & graph, istream & pepProph)
{
  BasicBigraph bb;
  bb.readFromMCMC(graph, pepProph);

  //cout << "Separating cleverly ;)" << endl;

  numberClones = 0;
  severedProteins = Array<string>();

#ifndef NOSEPARATE
#ifndef NOPRUNE
  Array<BasicBigraph> subBasic = iterativePartitionSubgraphs(bb, 0.0);
#else
  Array<BasicBigraph> subBasic = iterativePartitionSubgraphs(bb, -1.0);  
#endif
#else
  Array<BasicBigraph> subBasic = Array<BasicBigraph>(1, bb);
#endif

  subgraphs = Array<BasicGroupBigraph>(subBasic.size());

  for (int k=0; k<subBasic.size(); k++)
    {
      subgraphs[k] = BasicGroupBigraph(subBasic[k]);
    }

  initialize();
}

double GroupPowerBigraph::getLogNumberStates() const
{
  double total = 0;
  for (int k=0; k<subgraphs.size(); k++)
    {
      total = Numerical::logAdd(total, subgraphs[k].logNumberOfConfigurations());
    }

  // add one because each peptide needs to be estimated once
  return total + 1;
}

Array<BasicBigraph> GroupPowerBigraph::iterativePartitionSubgraphs(BasicBigraph & bb, double newPeptideThreshold )
{
  //  cerr << "Iter partition... @ threshold = " << newPeptideThreshold << endl;
  bb.PeptideThreshold = newPeptideThreshold;

  bb.prune();

  severedProteins.append( bb.severedProteins );
  numberClones += bb.numberClones;

  Array<BasicBigraph> preResult = bb.partitionSections();
  Array<BasicBigraph> result;

  //  cerr << "Using threshold " << newPeptideThreshold << " split into " << preResult.size() << endl;
  for (int k=0; k<preResult.size(); k++)
    {
      BasicGroupBigraph bgb = BasicGroupBigraph( preResult[k] );
      double logNumConfig = bgb.logNumberOfConfigurations();
      if ( logNumConfig > LOG_MAX_ALLOWED_CONFIGURATIONS && log2(bgb.PSMsToProteins.size())+log2(bgb.originalN[0].size+1) <= LOG_MAX_ALLOWED_CONFIGURATIONS )
	{
	  // the graph can, theoretically, become pruned to the
	  // desired efficiency, but has not been yet

	  //	  double smallest = Vector(preResult[k].PSMsToProteins.weights).min();
	  double newThresh = 1.25*(newPeptideThreshold + 1e-6);

	  Array<BasicBigraph> completelyFragmented = iterativePartitionSubgraphs(preResult[k], newThresh);

	  result.append( completelyFragmented );
	}
      else if (logNumConfig > LOG_MAX_ALLOWED_CONFIGURATIONS)
	{
	  // the graph cannot become pruned to the desired efficiency;
	  // prune as much as possible
	  double largest = Vector(preResult[k].PSMsToProteins.weights).max();

	  Array<BasicBigraph> completelyFragmented = iterativePartitionSubgraphs(preResult[k], largest);

	  result.append( completelyFragmented );
	}
      else
	{
	  // the graph is already pruned to the desired degree
	  result.add( preResult[k] );
	}
    }

  return result;
}

#ifndef NOSEPARATE
void GroupPowerBigraph::read(istream & is)
{
    // cout << "Reading GroupPowerBigraph" << endl;
    
  //  cout << "\tReading BasicBigraph" << endl;

  BasicBigraph bb;
  is >> bb;

  numberClones = 0;
  severedProteins = Array<string>();

  //  cout << "Starting partition..." << endl;

  #ifndef NOPRUNE
  Array<BasicBigraph> subBasic = iterativePartitionSubgraphs(bb, PeptideThreshold);
  #else
  Array<BasicBigraph> subBasic = iterativePartitionSubgraphs(bb, -1);
  #endif
  
  //  cout << "Finished partitioning into " << subBasic.size() << " partitions" << endl;

  // useful for debugging the iterative partition for small graphs
  //  for (int cc=0; cc<subBasic.size(); cc++)
  //    {
  //subBasic[cc].displayDotty("Test");
  //    }

  //  cout << "Number of clones was " << numberClones << endl << endl;

  //  cout << "\tPartition complete" << endl;

  //  cout << "Total number of proteins is " << bb.proteinsToPSMs.size() << endl;
  //  cout << "Total number of PSMs is " << bb.PSMsToProteins.size() << endl;

  //  exit(0);

  subgraphs = Array<BasicGroupBigraph>(subBasic.size());

  for (int k=0; k<subBasic.size(); k++)
    {
      subgraphs[k] = BasicGroupBigraph(subBasic[k]);
    }

  initialize();
  //  cout << "Initialized" << endl;
}

#else

 // for brute force: put into one subgraph
 void GroupPowerBigraph::read(istream & is)
 {
 BasicBigraph bb;
 is >> bb;

 bb.prune();

 numberClones = bb.numberClones;

 subgraphs = Array<BasicGroupBigraph>(1, BasicGroupBigraph(bb) );
 initialize();
 }
#endif

void GroupPowerBigraph::initialize()
{
  getGroupProtNames();
  //  getSumLikelihoodOverAllAlphaBetaGivenD();

  //  cout << "Initialized" << endl;
}

/***
void GroupPowerBigraph::gridScan() const
{
  DetectabilityModel local = gm;

  Array<Array<double> > mat( local.maxRowIndex(), Array<double>( local.maxColIndex(), -1) );
  double bestProb = -1;
  Model bestModel;
  for (local.start(); local.inRange(); local.advance())
    {
      double prob = probabilityAlphaBetaGivenD(local);
      //      double prob = logLikelihoodAlphaBetaGivenD(local);

      if ( bestProb < prob )
	{
	  bestProb = prob;
	  bestModel = Model(local);

	  //	  cout << "\tFound new best " << bestProb << endl << bestModel << endl;
	}

      mat[ local.getRowIndex() ][ local.getColIndex() ] = prob;
      //      mat[ local.getRowIndex() ][ local.getColIndex() ] = unprunedProbabilityAlphaBetaGivenD(local);
    }

  cerr << "Best model: " << bestModel << endl;
  cerr << "\thad prob: " << bestProb << endl;

  Matrix(mat).displayMatrix();
}
***/

ostream & operator <<(ostream & os, pair<double,double> rhs)
{
  os << "("<< rhs.first << ", " << rhs.second << ")";
  return os;
}

void GroupPowerBigraph::outputPivdo(ostream & os) const
{
  for (int k=0; k<subgraphs.size(); k++)
    PivdoSplitter(subgraphs[k]).outputPivdo(os);
}

Array<string> GroupPowerBigraph::peptideNames() const
{
  Array<string> pepNames;

  for (int k=0; k<subgraphs.size(); k++)
    {
      pepNames.append( subgraphs[k].PSMsToProteins.names );
    }

  return pepNames;
}

void GroupPowerBigraph::loadPeptideDetectability(char*fname)
{
  //  cout << "reading pep detect." << endl;

  cerr << "Reading PepFtrReader" << endl;
  pepFtrRdr.read(fname);
  //  cerr << "Reordering PepFtrReader" << endl;
  //  pepFtrRdr.reorder(peptideNames());
  //  cout << "pep detect has been read." << endl;
}

void GroupPowerBigraph::showDetectabilities()
{
  for (int k=0; k<subgraphs.size(); k++)
    {
      for (int j=0; j<subgraphs[k].proteinsToPSMs.size(); j++)
	{
	  cout << subgraphs[k].proteinsToPSMs.names[j] << "\t";
	  const Set & s = subgraphs[k].proteinsToPSMs.associations[j];
	  
	  cout << subgraphs[k].PSMsToProteins.names[ s ] << "\t";
	  cout << subgraphs[k].PSMsToProteins.weights[ s ] << endl;
	  /***
	  for (Set::Iterator iter = s.begin(); iter != s.end(); iter++)
	    {
	      cout << "\t" << subgraphs[k].PSMsToProteins.names[*iter] << "  \t alpha = " << subgraphModels[k].alphaArray[*iter] << endl;
	    }
	  ***/
	}
    }
}

void GroupPowerBigraph::loadKnownDetectabilities(char * fname)
{
  cerr << "Sorting subgraph peptide names for n log(n) lookup" << endl;
  Array<Array<string> > sortedSubgraphPeps(subgraphs.size());
  Array<Array<int> > sortedSubgraphIndices(subgraphs.size());

  for (int k=0; k<subgraphs.size(); k++)
    {
      sortedSubgraphPeps[k] = subgraphs[k].PSMsToProteins.names;
      sortedSubgraphIndices[k] = sortedSubgraphPeps[k].sortA();
    }

  //  cout << sortedSubgraphPeps << endl;
  //  cout << sortedSubgraphIndices << endl << endl;

  cerr << "Loading known detectabilities" << endl;
  ifstream fin(fname);
  
  string pep;
  double prob;
  while (fin >> pep >> prob)
    {
      changePeptideDetectability(sortedSubgraphPeps, sortedSubgraphIndices, pep, prob);
    }
}

void GroupPowerBigraph::loadKnownDetectabilitiesFast(char * fname)
{
  cerr << "Loading known detectabilities" << endl;
  ifstream fin(fname);
  
  Array<string> allPeps;
  Array<double> allProbs;

  string pep;
  double prob;
  while (fin >> pep >> prob)
    {
      allPeps.add(pep);
      allProbs.add(prob);
    }

  cerr << "Sorting all peptide names for n log(n) lookup" << endl;

  Array<string> allPepsSorted = allPeps;
  Array<int> sortingMap = allPepsSorted.sortA();
  Array<double> allProbsSorted = allProbs[sortingMap];

  cerr << "Setting values" << endl;

  //  cerr << allPepsSorted[0] << endl;

  double alphaInitial = subgraphModels[0].alphaArray[0];
  cerr << "Using alphaInitial " << alphaInitial << endl;
  for (int k=0; k<subgraphs.size(); k++)
    {
      for (int j=0; j<subgraphs[k].PSMsToProteins.size(); j++)
	{
	  const string & pepName = subgraphs[k].PSMsToProteins.names[j];
	  double & alpha = subgraphModels[k].alphaArray[j];

	  int ind = allPepsSorted.sortedFind(pepName);
	  if ( ind != -1 )
	    {
	      //	      alpha = clip(0.08);
	      alpha = clip( allProbsSorted[ind] * alphaInitial );
	      //	      cerr << pepName << "\t" << weight << endl;
	    }
	  else
	    {
	      // when the peptide isn't found in the detectability
	      // file, let it be very undetectable
	      alpha = 0.0001;
	    }
	}
    }
}

void GroupPowerBigraph::changePeptideDetectability(const Array<Array<string> > & sortedSubgraphPeps, const Array<Array<int> > & sortedSubgraphIndices, const string & pepStr, double prob)
{
  //  cout << "Index of peptide " << pepStr << endl;
  for (int k=0; k<subgraphs.size(); k++)
    {
      int rearrangedIndex = sortedSubgraphPeps[k].sortedFind(pepStr);

      //      cout << "\tsortedFind res for " << sortedSubgraphPeps[k] << " is " << sortedSubgraphPeps[k].sortedFind(pepStr) << endl;
      if ( rearrangedIndex != -1 )
	{
	  int pepIndex = sortedSubgraphIndices[k][ rearrangedIndex ];
	  //	  cout << "\t(@subgraph " << k << ") is " << pepIndex << endl;
	  subgraphModels[k].alphaArray[pepIndex] = clip(prob);

	  //	  cerr << pepStr << " " << clip(prob);
	}
    }
}

void GroupPowerBigraph::setAlphaBetaGamma(double alpha, double beta, double gamma)
{
  subgraphModels = Array<DetectabilityModel>(subgraphs.size());
  for(int k=0; k<subgraphModels.size(); k++)
    {
      subgraphModels[k] = DetectabilityModel(Array<double>(subgraphs[k].PSMsToProteins.size(), alpha), beta, gamma);
    }  
}

double GroupPowerBigraph::clip(double x) const
{
  // for numerical stability... grr
  return min(.999, max(1e-3, x) );
}

