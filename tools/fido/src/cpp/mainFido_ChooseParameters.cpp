
#include <climits>

#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <set>

#include <unistd.h>

#include "GroupPowerBigraph.h"
#include "InputFileShrinker.h"


#include <time.h>
#include <assert.h>

using namespace std;

double antiderivativeAt(double m, double b, double xVal)
{
  return m*xVal*xVal/2.0 + b*xVal;
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

double squaredArea(double x1, double y1, double x2, double y2, double max_x)
{
  if ( x2 < x1 )
    return 0.0;
  assert(x2 >= x1);
  double m = (y2-y1)/(x2-x1);
  double b = y1-m*x1;

  double result = squareAntiderivativeAt(m, b, min(max_x, x2) ) - squareAntiderivativeAt(m, b, x1);
  //  if (result < 0.0)
  //  std::cout << "sq_area " << result << std::endl << m << " " << b << " " << x1 << " " << x2 << std::endl;;
  return result;
}

double area(double x1, double y1, double x2, double y2, double max_x)
{
  //  std::cout << x1 << " " << x2 << std::endl;
  assert(x2 >= x1);
  double m = (y2-y1)/(x2-x1);
  double b = y1-m*x1;

  double result = antiderivativeAt(m, b, min(max_x, x2) ) - antiderivativeAt(m, b, x1);
  if (result < 0.0)
    std::cout << "area " << result << std::endl << m << " " << b << " " << x1 << " " << x2 << std::endl;;
  return result;
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
      cerr << "Warning: There are not enough false positives; needed " << N << " and was only given " << fpArray.back() << " (will proceed using largest available value)" << endl << endl;
      //      exit(1);
      N = fpArray.back();
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
  double totalFDR = 0.0, estFDR = 0.0, empFDR;

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

	  //	  std::cout << "est:" << estFDR << " " << totalFDR << " " << (1-prob) << " " << (fpChange+tpChange) << (fpCount + tpCount) << std::endl;
	  //	  std::cout << prob << " " << 1-prob << " " << estFDR << " " << empFDR << std::endl;
	  
	  if ( estFDR > estFDR_array[estFDR_array.size()-1] )
	    {
	      estFDR_array.add(estFDR);
	      empFDR_array.add(empFDR);
	    }
	}
      totalFDR += (1-prob) * (fpChange + tpChange);
      fpCount += fpChange;
      tpCount += tpChange;
      estFDR = totalFDR / (fpCount + tpCount);

      empFDR = double(fpCount) / (fpCount + tpCount);
      lastProb = prob;
    }
  
  return pair<Array<double>, Array<double> >(estFDR_array, empFDR_array);
}

double getFDR_divergence(const Array<double> estFDR, const Array<double> empFDR, double THRESH)
{
  Vector diff = Vector(estFDR) - Vector(empFDR);
  //  std::cout << "est " << estFDR << std::endl;
  //  std::cout << "emp " << empFDR << std::endl;

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

      // use estimated FDR as the x-axis (it is guaranteed to be
      // sorted, while empirical FDR is not)
      tot += squaredArea(estFDR[k], diff[k], estFDR[k+1], diff[k+1], estFDR[k+1]);
    }

  //  std::cout << "tot: " << tot << std::endl;
  double xRange = min(THRESH, estFDR[k]) - estFDR[0];
  //  std::cout << "xr: " << xRange << std::endl;

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

void print_usage()
{		
		cerr << "usage: FidoChooseParameters [-p] [-a] [-g] [-c <n>] <graph file> <target decoy file>" << endl;
        cerr << "       FidoChooseParameters [-p] [-a] [-g] [-c <n>] <graph file> <target decoy file> <l1>" << endl;
        cerr << "       FidoChooseParameters [-p] [-a] [-g] [-c <n>] <graph file> <target decoy file> <l1> <l2>" << endl;
        
		cerr << "where" << endl;
		cerr << "          l1 is the log2 of maximum number of subgraph connected states, " << endl;
		cerr << "          l2 is the log2 for the main calculation, and l1 is only used , " << endl;
		cerr << "             for precalculation," << endl;
		cerr << "          option -p omits cleaning the peptide names," << endl;
		cerr << "          option -a uses all PSM matches instead the best one, and" << endl;
		cerr << "          option -g uses protein group level inference." << endl;
		cerr << "          option -c sets start parameter's accurary level (1-3)" << endl;
		cerr << "                                        1 = best    / slower" << endl;
		cerr << "                                        2 = relaxed / faster" << endl;
		cerr << "                                        3 = sloppy  / very fast" << endl;
}

// global vars used to store command line options during runtime
extern bool gbDoPeptideNameClearing;
extern bool gbUseAllPepMatches;
extern bool gbUseProteinGroupLevelInference;
int  giParameterAccuracySetting;
string  gsrReducedParameterFilename;

void set_default_cmdline_opts()
{
	giParameterAccuracySetting = 1; // use full file
	gbDoPeptideNameClearing = true;
	gbUseAllPepMatches = false;
	gbUseProteinGroupLevelInference = false;
}



int main(int argc, char**argv)
{
	int option_char;
	char *s;
	
	// set default command line options
	set_default_cmdline_opts();
	
	while ((option_char = getopt (argc, argv, "pagc:" )) != EOF)
    switch (option_char)
      {  
		case 'p': cerr << "omitting peptide sequence cleaning" << endl; 
			gbDoPeptideNameClearing = false;
			break;

		case 'a': cerr << "use all peptide matches" << endl; 
			gbUseAllPepMatches = true;
			break;
		case 'g': cerr << "use protein group level inference" << endl; 
			gbUseProteinGroupLevelInference = true;
			break;
		case 'c' : 
				s = optarg;
				if(s == NULL){
					print_usage();
					return 1;
				}
				
				giParameterAccuracySetting = atoi(optarg);
				
				cerr << "attempting to use parameter accuracy settting: " << giParameterAccuracySetting << endl;
        
		 default: 
			 break;
      }

try{

// 2: only filenames given
// 3: additional log2 parameter given
// 4: two log2 parameters given
  if ( argc - (optind) == 2 || argc - (optind) == 3 || argc - (optind) == 4 )
    {	
	  //cout << "argc = " << argc << endl;
	  //cout << "optind = " << optind << endl;
	  //cout << "input file argv[0 + optind] = " << argv[0 + optind] << endl;
	  //cout << "decoy file argv[1 + optind] = " << argv[1 + optind] << endl;
	  //exit(0);
		
      //srand(time(NULL));
      cout.precision(10);

      char*tdecoy_fname = argv[1+(optind)];
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

      if ( argc - (optind) >= 3 )
		GroupPowerBigraph::LOG_MAX_ALLOWED_CONFIGURATIONS = atof(argv[2+(optind)]);

      //      cerr << "\tLoading..." << endl;
      double gamma, alpha, beta;
      gamma = alpha = beta = 0.01;

		
		
		// targetSize is the size of the dataset that is used to estimate $\alpha$, $\beta$ and $\gamma$.
	unsigned int targetSize;
		switch(giParameterAccuracySetting){
			case 1: targetSize = INT_MAX;   
				break;
			case 2: targetSize = 2500;
				break;
			case 3: targetSize = 100;
				break;
			case 42: targetSize = 30000;
				break;
			default: targetSize = INT_MAX;   
				print_usage();
				exit(1);
				break;
		}



	GroupPowerBigraph *gpbPrecalc;
	//GroupPowerBigraph gpbfull( argv[0+(optind)], gamma, alpha, beta);
	
		

		// todo: replace this with: if(targetSize < gpfull.NumberOfPeptides() 
		if(targetSize < INT_MAX){
			gsrReducedParameterFilename = (string(argv[0+(optind)])+ string("_shrink"));
		InputFileShrinker shrinker;
			shrinker.shrink(argv[0+(optind)], gsrReducedParameterFilename.c_str(), targetSize );
			gpbPrecalc = new GroupPowerBigraph( (char *)gsrReducedParameterFilename.c_str(), gamma, alpha, beta);	
		}
		else
		{
			gpbPrecalc = new GroupPowerBigraph( argv[0+(optind)], gamma, alpha, beta);	
		}

      //      cerr << "\tMarginalizing..." << endl;

      double gamma_best, alpha_best, beta_best;
      gamma_best = 0.5;
      alpha_best = 0.1;
      beta_best = 0.05;
      double best_objective = -100000000;

      double gamma_search[] = {0.1, 0.3, 0.5, 0.7, 0.9};
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
		  
		  
		  (*gpbPrecalc).setAlphaBetaGamma(alpha, beta, gamma);
		  (*gpbPrecalc).getProteinProbs();
		  pair<Array<Array<string> >, Array<double> > prot_groups_and_probs = (*gpbPrecalc).getDescendingProteinsAndWeights();
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
      
      // now for the main calculation 
    GroupPowerBigraph *gpbMaincalc;
	  if( (giParameterAccuracySetting == 1) && (argc - (optind) != 4) ){
		  gpbMaincalc = gpbPrecalc;
	  }else{
		  delete gpbPrecalc;
		  
		  if(argc - (optind) == 4){
			  GroupPowerBigraph::LOG_MAX_ALLOWED_CONFIGURATIONS = atof(argv[3+(optind)]);
		  }
		  gpbMaincalc = new GroupPowerBigraph( argv[0+(optind)], gamma, alpha, beta);
	  }
      
      (*gpbMaincalc).setAlphaBetaGamma(alpha_best, beta_best, gamma_best);
      (*gpbMaincalc).getProteinProbs();
      (*gpbMaincalc).printProteinWeights();
      
	  delete gpbMaincalc;      
    }
  else
    {
		print_usage();
		return 1;
    }
  }catch(exception &e)
  {
	  cerr << "caught an exception :" << e.what() << endl;
	  return -1;
  }
  catch(...){
	  cerr << "something strange just happened" << endl;
	  return -1;
  }

  return 0;
}

