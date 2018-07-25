#include "Proteomics.h"

int Proteomics::maxLength, Proteomics::minLength;

Array<string> Proteomics::trypsinDigest(const string & str)
{
  Array<string> result;

  int n = str.length();
  int loc = 0;
  while ( loc < n )
    {
      bool cleavageFound = false;

      for (int k=loc+1; k<n-1; k++)
	{
	  if ( ( str[k] == 'R' || str[k] == 'K' ) && str[k+1] != 'P' )
	    {
	      result.add( str.substr(loc, k-loc+1) );
	      loc = k+1;
	      cleavageFound = true;

	      break;
	    }
	}

      if ( ! cleavageFound )
	{
	  result.add( str.substr(loc, n-loc) );
	  break;
	}
    }

  return result;
}

bool Proteomics::observablePeptide(const string & str)
{
  return (int)str.length() <= maxLength && (int)str.length() >= minLength;
}

void Proteomics::setMinMaxLength(char * pivdoFname, double minScore)
{
  ifstream fin(pivdoFname);
  BasicBigraph bb;
  fin >> bb;

  const Array<string> & peptides = bb.PSMsToProteins.names;
  const Array<double> & probs = bb.PSMsToProteins.weights;

  Set goodProbsSet = Vector(probs) >= minScore;

  Array<double> goodProbs = probs[ goodProbsSet ];
  Array<int> goodPepLengths;
  for (int k=0; k<goodProbsSet.size(); k++)
    {
      goodPepLengths.add( peptides[ goodProbsSet[k] ].length() );
    }

  goodPepLengths.sort();

  // take only between the 5th and 95th percentile
  //  cout << goodPepLengths << endl;

  maxLength = goodPepLengths[ int(goodPepLengths.size() * .05) ];
  minLength = goodPepLengths[ int(goodPepLengths.size() * .95) ];

  cerr << "Min len = " << minLength << endl;
  cerr << "Max len = " << maxLength << endl;
}

