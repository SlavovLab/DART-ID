#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include "Vector.h"
#include <set>
//#include "Hash_Set.h"



using namespace std;

#include <ext/hash_set>

using namespace std;

namespace __gnu_cxx
{
  template<> struct hash< std::string >
  {
    size_t operator()( const std::string & x ) const
    {
      return hash< const char* >()( x.c_str() );
    }
  };
}


int matchCount( const __gnu_cxx::hash_set<string> & positiveNames, const Array<string> & atThreshold )
{
  int count = 0;
  
  for (int k=0; k<atThreshold.size(); k++)
    {
      if ( positiveNames.count( atThreshold[k] ) > 0 )
	count++;
    }

  return count;
}

Array<string> matches( const __gnu_cxx::hash_set<string> & positiveNames, const Array<string> & atThreshold )
{
  Array<string> result;
  for (int k=0; k<atThreshold.size(); k++)
    {
      if ( positiveNames.count( atThreshold[k] ) > 0 )
	result.add( atThreshold[k] );
    }
  return result;
}

int main(int argc, char**argv)
{
  if ( argc == 2 )
    {
      ifstream fin(argv[1]);

      Array<string> truePositiveNames, falsePositiveNames;
      fin >> truePositiveNames;
      fin >> falsePositiveNames;

      __gnu_cxx::hash_set<string> truePosSet(truePositiveNames.size()), falsePosSet(falsePositiveNames.size());

      int k;
      for (k=0; k<truePositiveNames.size(); k++)
	{
	  truePosSet.insert( truePositiveNames[k] );
	  //	  cout << "\tinsert " << truePositiveNames[k];
	}
      for (k=0; k<falsePositiveNames.size(); k++)
	{
	  falsePosSet.insert( falsePositiveNames[k] );
	  //	  cout << "\tinsert " << falsePositiveNames[k];
	}

      //      cout << "Found TP, FP names: " << endl;
      //            cout << truePositiveNames << endl;
      //      cout << falsePositiveNames << endl << endl;

      Array<string> protsAtThreshold;
      string line;

      double prob, lastProb=-1;

      Array<int> fps;
      Array<Array<string> > tps;

      int fpCount = 0, tpCount = 0;

      int numScored = 0;
      Array<string> observedProteins;

      bool scheduledUpdate = false;

      while ( cin >> prob && getline(cin, line) )
	{
	  //	  cerr << "prob = " << prob << ", \tline = " << line << endl;

	  istringstream ist(line);
	  ist >> protsAtThreshold;

	  //	  cerr << "\tRead " << protsAtThreshold << endl;

	  numScored += protsAtThreshold.size();

	  int fpChange = matchCount(falsePosSet, protsAtThreshold);
	  int tpChange = matchCount(truePosSet, protsAtThreshold);

	  protsAtThreshold.sort();
	  protsAtThreshold = Array<string>(1, protsAtThreshold[0]);
	  if ( fpChange > 0 && tpChange > 0)
	    {
	      fpChange = tpChange = 0;
	      protsAtThreshold = Array<string>();
	    }
	  else if ( fpChange > 0 )
	    fpChange = 1;
	  else if ( tpChange > 0 )
	    tpChange = 1;

	  //	  cout << "This iter has fpChange, tpChange = " << fpChange << ", " << tpChange << endl;
	  //	  cout << "\tNumber fps, tps = " << fpChange << ", " << tpChange << endl;

	  if ( prob != lastProb && lastProb != -1 )
	    {
	      scheduledUpdate = true;
	      //	      cout << "prob has changed from " << lastProb << " to " << prob << endl;
	    }
	  
	  if ( scheduledUpdate )
	    {
	      //	      if ( fpChange > 0 )
		{
		  //		  cout << fpCount << " ";
		  //		  cout << matches(truePosSet, observedProteins) << endl;

		  //		  if ( fpCount+fpChange > 10 )
		  //		    {
		  //  		      cout << " " << matches(falsePosSet, observedProteins) << endl;
		  //		      exit(0);
		  //		    }

		  fps.add( fpCount );
		  //		  tps.add( matches(truePosSet, observedProteins) );
		  tps.add( matches(truePosSet, observedProteins) );

		  //		  cout << fpCount << " " << matchCount(truePosSet, observedProteins) << endl;

		  //		  cout << endl;
		  //		  cout << "Adding cumulative totals " << fpCount << ", " << tpCount << endl;
		  scheduledUpdate = false;

		  if ( fpCount > 100 )
		    break;

		  if ( tpChange > 0 )
		    {
		      //		      		      cerr << "Adding new TPs at threshold " << prob << endl;
		      //		      		      cerr << "\tto make data point with FP, TP = " << fpCount+fpChange << ", " << tpCount+tpChange << endl;
		      //		      		      cerr << "\tTPs added were " << matches(truePosSet, protsAtThreshold ) << endl;
		    }
		}
	    }

	  fpCount += fpChange;
	  tpCount += tpChange;


	  /***
	  if ( fpCount == 10 )
	    {
	      Array<string> tpsAt10 = matches(truePosSet, observedProteins);
	      cout << "AT10: " << tpsAt10 << endl;
	      return 0;
	    }
	  ***/


	  lastProb = prob;
	  observedProteins.append( protsAtThreshold );
	}

      //      cerr << "Fin loop" << endl;

      lastProb = prob;

      //      if ( scheduledUpdate )
	{
	  //	  fps.add( fpCount );
	  //	  tps.add( tpCount );
	}

	//		cerr << "Final counts are " << fpCount << ", " << tpCount << endl;
	//		cerr << "Overall: " << matchCount(falsePosSet, observedProteins) << ", " << matchCount(truePosSet, observedProteins) << endl;
	//		cerr << "Total number of proteins scored is " << numScored << endl;

	//	cout << "end" << endl << endl;
	//	cout << falsePosSet.size() << " " << truePositiveNames << endl;

	fps.add( falsePosSet.size() );
	tps.add( truePositiveNames );
	
	//	cout << "------------------" << endl;
	
	/***
	for (int u=0; u<fps.size(); u++)
	  {
	    cout << fps[u] << " " << tps[u] << endl;
	  }

	cout << endl;
	cout << "------------------" << endl;
	***/

	cout << fps << " " << tps << endl;
    }
  else
    {
      cerr << "usage: <True Positives and Overall Array File>" << endl << endl;
    }
  return 0;
}


