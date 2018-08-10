#ifndef _GroupPowerBigraph_H
#define _GroupPowerBigraph_H

#include "ProteinIdentifier.h"
#include "BasicGroupBigraph.h"
#include "PivdoSplitter.h"
#include "PeptideFeatureReader.h"

using namespace std;


class GroupPowerBigraph : public ProteinIdentifier
{
public:
  static double LOG_MAX_ALLOWED_CONFIGURATIONS;

 GroupPowerBigraph(char * fname, double gamma, double alpha, double beta) :
  ProteinIdentifier(),
  pepFtrRdr(1)
      {
	//	cerr << "Loading without an rbfGamma and using a pruning threshold of " << LOG_MAX_ALLOWED_CONFIGURATIONS << endl;

	ifstream fin(fname);
	read(fin);
	
	setAlphaBetaGamma(alpha, beta, gamma);
      }

 GroupPowerBigraph(char * fname, char * detFname, double rbfGamma, double gamma, double alphaInitial, double beta):
  ProteinIdentifier(),
  pepFtrRdr(rbfGamma)
      {
	//	cerr << "Loading with rbfGamma = " << rbfGamma << " and using a pruning threshold of " << LOG_MAX_ALLOWED_CONFIGURATIONS << endl;

	ifstream fin(fname);
	read(fin);
	
	cerr << "Loading detect" << endl;
	
	loadPeptideDetectability(detFname);

	cerr << "Making subgraph models" << endl;

	subgraphModels = Array<DetectabilityModel>(subgraphs.size());

	for(int k=0; k<subgraphModels.size(); k++)
	  {
	    subgraphModels[k] = DetectabilityModel(Array<double>(subgraphs[k].PSMsToProteins.size(), alphaInitial), beta, gamma);
	    //	    cerr << "Constructed with model " << subgraphModels[k].alphaArray << " " << subgraphModels[k] << endl;
	  }
      }

  void setAlphaBetaGamma(double alpha, double beta, double gamma);

  double clip(double) const;

  void printGraphs();

  void showDetectabilities();
  void printPeptideDetectabilities();
  void loadKnownDetectabilitiesFast(char * fname);
  void loadKnownDetectabilities(char * fname);
  void changePeptideDetectability(const Array<Array<string> > & sortedSubgraphPeps, const Array<Array<int> > & sortedSubgraphIndices, const string & pepStr, double prob);

  Array<double> proteinProbs();
  Array<double> MAPProteins();
  void printProteinWeights() const;

  void getProteinProbs();
  void getMAPProteins();

  Array<string> peptideNames() const;
  void scanAlphaBeta();

  Array<BasicBigraph> iterativePartitionSubgraphs(BasicBigraph & bb, double newPeptideThreshold );

  double getLogNumberStates() const;

  double logLikelihoodAlphaBetaGivenD() const;
  /***
  Array<double> proteinProbsOverAllAlphaBeta();
  void getProteinProbsOverAllAlphaBeta();
  double probabilityAlphaBetaGivenD();

  double likelihoodAlphaBetaGivenD( const DetectabilityModel & myGM ) const;

  // only pass this so that it can be cached easily 
  // even though you are passing a member variable
  double sumLogLikelihoodOverAllAlphaBeta( const DetectabilityModel & myGM ) const;

  // only pass this so that it can be cached easily 
  // even though you are passing a member variable

  double probabilityAlphaBetaGivenD( const DetectabilityModel & myGM ) const;

  void gridScan() const;
  **/

  void readFromMCMC(istream & graph, istream & pepProph);

  void outputPivdo(ostream & os) const;

  void loadPeptideDetectability(char*fname);

  void refreshModels();
  double meanAlpha(const Array<double> & probXAndEGivenD, const Array<double> & probXGivenD);

  Array<Array<string> > groupProtNames;

  pair<Array<Array<string> >, Array<double> > getDescendingProteinsAndWeights() const;

protected:

  Array<double> probabilityXAndEGivenD();
  Array<double> probabilityXGivenD();

  Array<string> severedProteins;

  int numberClones;
  void initialize();

  void getProteinProbsGivenAlphaBeta();

  void getGroupProtNames();

  Array<double> probabilityR;

  void read(istream & is);

  Array<BasicGroupBigraph> subgraphs;
  Array<DetectabilityModel> subgraphModels;

  // cached functors
  //  LastCachedMemberFunction<GroupPowerBigraph, double, DetectabilityModel> sumLogLikelihoodOverAllAlphaBetaCachedFunctor;

  Numerical zeroChecker;
  PeptideFeatureReader pepFtrRdr;
};


ostream & operator <<(ostream & os, pair<double,double> rhs);


#endif

