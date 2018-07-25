#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <set>

#include "GroupPowerBigraph.h"

using namespace std;

double antiderivativeAt(double m, double b, double xVal)
{
  return m*xVal*xVal/2.0 + b*xVal;
}

double area(double x1, double y1, double x2, double y2, double max_x)
{
  double m = (y2-y1)/(x2-x1);
  double b = y1-m*x1;

  return antiderivativeAt(m, b, min(max_x, x2) ) - antiderivativeAt(m, b, x1);
}

double squareAntiderivativeAt(double m, double b, double xVal)
{
  // turn into ux^2+vx+t
  double u = m*m;
  double v = 2*m*b;
  double t = b*b;

  //  cout << "\t\tsquareAntiderivativeAt " << xVal << " is " << u*xVal*xVal*xVal/3.0 + v*xVal*xVal/2.0 + t*xVal << endl;
  //  cout << "u, v, t = " << u << " " << v << " " << t << endl;

  return u*xVal*xVal*xVal/3.0 + v*xVal*xVal/2.0 + t*xVal;
}

int matchCount( const set<string> & positiveNames, const Array<string> & atThreshold )
{
  int count = 0;
  
  for (int k=0; k<atThreshold.size(); k++)
    {
      if ( positiveNames.count( atThreshold[k] ) > 0 )
	count++;
    }

  return count;
}

double getROC_N(const Array<int> & fpArray, const Array<int> & tpArray, int N)
{
  double rocN = 0.0;

  if ( fpArray.back() < N )
    {
      cerr << "There are not enough false positives; needed " << N << " and was only given " << fpArray.back() << endl << endl;
      exit(1);
    }

  for (int k=0; k<fpArray.size()-1; k++)
    {
      // find segments where the fp value changes
	  
      if ( fpArray[k] >= N )
	break;

      if ( fpArray[k] != fpArray[k+1] )
	{
	  // this line segment is a function
	      
	  double currentArea = area(fpArray[k], tpArray[k], fpArray[k+1], tpArray[k+1], N);
	  rocN += currentArea;
	}
    }
  return rocN / (N * tpArray.back());
}

pair<Array<double>, Array<double> > getEstimated_and_Empirical_FDR(Array<Array<string> > names, Array<double> probabilities, const set<string> & falsePosSet, const set<string> & truePosSet)
{
  Array<double> estFDR_array, empFDR_array;
  estFDR_array.add(0);
  empFDR_array.add(0);
  
  int fpCount, tpCount;
  fpCount = tpCount = 0;
  double totalFDR = 0.0, estFDR = 0.0;

  bool scheduledUpdate = false;
  
  double lastProb = -1.0;
  for (int k=0; k<names.size(); k++)
    {
      double prob = probabilities[k];
      int fpChange = matchCount(falsePosSet, names[k]);
      int tpChange = matchCount(truePosSet, names[k]);
      
      // for different style of grading, counting groups as a
      // single protein and throwing away any groups that include
      // TPs and FPs
      if ( tpChange > 0 && fpChange > 0 )
	tpChange = fpChange = 0;
      
      if ( tpChange > 0 )
	tpChange = 1;
      if ( fpChange > 0 )
	fpChange = 1;
      
      if ( prob != lastProb && lastProb != -1 )
	{
	  scheduledUpdate = true;
	}

      if ( scheduledUpdate )
	{
	  scheduledUpdate = false;

	  totalFDR += (1-prob) * (fpChange + tpChange);
	  estFDR = totalFDR / (fpCount + tpCount);
	  
	  estFDR_array.add(estFDR);
	  empFDR_array.add(totalFDR);
	}

      fpCount += fpChange;
      tpCount += tpChange;

      lastProb = prob;
    }
  
  return pair<Array<double>, Array<double> >(estFDR_array, empFDR_array);
}

double getFDR_divergence(const Array<double> estFDR, const Array<double> empFDR, double THRESH)
{
  Vector diff = Vector(estFDR) - Vector(empFDR);

  double tot = 0.0;

  int k;
  for (k=0; k<diff.size()-1; k++)
    {
      // stop if no part of the estFDR is < threshold
      if ( estFDR[k] >= THRESH )
	{
	  if ( k == 0 )
	    tot = 1.0 / 0.0;

	  break;
	}

      tot += area(estFDR[k], diff[k], estFDR[k+1], diff[k+1], estFDR[k+1]);
    }

  double xRange = min(THRESH, estFDR[k]) - estFDR[0];

  if ( isinf(tot) )
    return tot;

  return tot / xRange;
}

pair<Array<int>, Array<int> > getROC(Array<Array<string> > names, Array<double> probabilities, const set<string> & falsePosSet, const set<string> & truePosSet)
{
  Array<int> fps, tps;
  fps.add(0);
  tps.add(0);
  
  int fpCount, tpCount;
  fpCount = tpCount = 0;

  bool scheduledUpdate = false;
  
  double lastProb = -1.0;
  for (int k=0; k<names.size(); k++)
    {
      double prob = probabilities[k];
      int fpChange = matchCount(falsePosSet, names[k]);
      int tpChange = matchCount(truePosSet, names[k]);
      
      // for different style of grading, counting groups as a
      // single protein and throwing away any groups that include
      // TPs and FPs
      if ( tpChange > 0 && fpChange > 0 )
	tpChange = fpChange = 0;
      
      if ( tpChange > 0 )
	tpChange = 1;
      if ( fpChange > 0 )
	fpChange = 1;
      
      if ( prob != lastProb && lastProb != -1 )
	{
	  scheduledUpdate = true;
	}

      if ( scheduledUpdate )
	{
	  fps.add( fpCount );
	  tps.add( tpCount );
	  scheduledUpdate = false;

	  //	  cout << fpCount << " " << tpCount << endl;

	  //	  totalFDR += (1-prob) * (fpChange + tpChange);
	  //	  estFDR = totalFDR / (fpCount + tpCount);
	}

      fpCount += fpChange;
      tpCount += tpChange;

      lastProb = prob;
    }

  fps.add( fpCount );
  tps.add( tpCount );	  
  
  fps.add( falsePosSet.size() );
  tps.add( truePosSet.size() );
  
  return pair<Array<int>, Array<int> >(fps, tps);
}

int main(int argc, char**argv)
{
  if ( argc == 3 || argc == 4 )
    {
      srand(time(NULL));
      cout.precision(8);

      char*tdecoy_fname = argv[2];
      ifstream fin(tdecoy_fname);
      Array<string> truePositiveNames, falsePositiveNames;
      fin >> truePositiveNames;
      fin >> falsePositiveNames;

      set<string> truePosSet, falsePosSet;

      int k;
      for (k=0; k<truePositiveNames.size(); k++)
	{
	  truePosSet.insert( truePositiveNames[k] );
	  //	  cout << "\tinsert " << truePositiveNames[k];
	}
      for (k=0; k<falsePositiveNames.size(); k++)
	{
	  falsePosSet.insert( falsePositiveNames[k] );
	}

      if ( argc == 4 )
	GroupPowerBigraph::LOG_MAX_ALLOWED_CONFIGURATIONS = atof(argv[3]);

      //      cerr << "\tLoading..." << endl;
      double gamma, alpha, beta;
      gamma = alpha = beta = 0.01;

      GroupPowerBigraph gpb( argv[1], gamma, alpha, beta);

      //      cerr << "\tMarginalizing..." << endl;

      double gamma_best, alpha_best, beta_best;
      gamma_best = alpha_best = beta_best = -1.0;
      double best_objective = -100000000;

      double gamma_search[] = {0.1, 0.5, 0.9};
      double alpha_search[] = {0.01, 0.04, 0.09, 0.16, 0.25, 0.36};
      double beta_search[] = {0.0, 0.01, 0.025, 0.05};
      for (unsigned int i=0; i<sizeof(gamma_search)/sizeof(double); i++)
	{
	  for (unsigned int j=0; j<sizeof(alpha_search)/sizeof(double); j++)
	    {
	      for (unsigned int k=0; k<sizeof(beta_search)/sizeof(double); k++)
		{
		  gamma = gamma_search[i];
		  alpha = alpha_search[j];
		  beta = beta_search[k];
		  gpb.setAlphaBetaGamma(alpha, beta, gamma);
		  gpb.getProteinProbs();
		  pair<Array<Array<string> >, Array<double> > prot_groups_and_probs = gpb.getDescendingProteinsAndWeights();
		  Array<Array<string> > prot_names = prot_groups_and_probs.first;
		  Array<double> prot_probs = prot_groups_and_probs.second;
		  pair<Array<int>, Array<int> > roc = getROC(prot_names, prot_probs, falsePosSet, truePosSet);
		  double roc50 = getROC_N(roc.first, roc.second, 50);
		  pair<Array<double>, Array<double> > est_and_emp_fdr = getEstimated_and_Empirical_FDR(prot_names, prot_probs, falsePosSet, truePosSet);
		  double fdr_mse = getFDR_divergence(est_and_emp_fdr.first, est_and_emp_fdr.second, 0.10);
		  double lambda = 0.15;
		  double current_objective = lambda * roc50 - (1-lambda) * fdr_mse;
		  if (current_objective > best_objective)
		    {
		      best_objective = current_objective;
		      gamma_best = gamma;
		      alpha_best = alpha;
		      beta_best = beta;
		    }
		  cerr << gamma << " " << alpha << " " << beta << " : " << roc50 << " " << fdr_mse << " " << current_objective << endl;
		}
	    }
	}

      cerr << "Using best gamma, alpha, beta = " << gamma_best << " " << alpha_best << " " << beta_best << endl;
      gpb.setAlphaBetaGamma(alpha_best, beta_best, gamma_best);
      gpb.getProteinProbs();
      gpb.printProteinWeights();
    }
  else
    {
      cerr << "usage: FidoChooseParameters <graph file> <target decoy file>" << endl;
      cerr << "       FidoChooseParameters <graph file> <target decoy file> <log2 of maximum number of subgraph connected states>" << endl;
    }
  return 0;
}

