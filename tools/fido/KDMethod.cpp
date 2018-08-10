#include "KDMethod.h"

void KDMethod::blockLinkSimilar(const Array<int> & from, const Array<int> & to)
{
  /***
  if ( from.size() * to.size() < 10000 )
    {
      cout << "Block link similar:" << endl;
      cout << "\t" << from << endl;
      cout << "x\t" << to << endl << endl;
    }
  ***/
  
  Array<int> cat = (Array<int>)from;
  cat.append( (Array<int>)to );

  Array<int> catFromBools = Array<int>(from.size(), 1);
  catFromBools.append( Array<int>(to.size(), 0) );

  // assumes that "to" does not have duplicates
  //  const Vector & centerVector = meanVector(to);

  // choose a random "from" point as the center
  const Vector & centerVector = *orderedFeatureVectorRefs[ from[ Random::inRange(0, from.size())] ];
  
  Array<double> prox(cat.size());
  int k;
  for (k=0; k<cat.size(); k++)
    prox[k] = diffNormSquared( centerVector, *orderedFeatureVectorRefs[cat[k]] );

  Array<int> sortedIndices = prox.sortA();
  catFromBools = catFromBools[ sortedIndices ];
  cat = cat[ sortedIndices ];

//   cout << prox << endl;
//   cout << cat << endl;
//   cout << catFromBools << endl;

  Array<int> beg = blockBeginnings(prox);

  cout << "blocks begin at " << beg << endl;
  cout << "epsilon is " << sparsify.epsilon << endl;
  cout << "number of prox is " << prox.size() << endl;

//   cout << "sortedIndices.size() == " << sortedIndices.size() << endl;
//   cout << "prox.size() == " << prox.size() << endl;
//   cout << "Split into " << beg.size() << " close blocks" << endl;

  if ( beg.size() == 1 )
    {
    // all must be close
      //      cerr << "Warning, this problem was split into only one collection of " << from.size() << "x" << to.size() << " points" << endl;
      linkSimilar(from, to, true);
    }
  else
    {
      for (k=0; k<beg.size()-1; k++)
	{
	  int secondLow = beg[k+1];
	  int firstHigh = secondLow - 1;
	  int firstLow = beg[k];

	  int secondHigh;
	  if ( k != beg.size() - 2 )
	    secondHigh = beg[k+2]-1;
	  else
	    secondHigh = prox.size()-1;

	  //      cout << "\tblock " << k << " contains " << cat[Set::FullSet(firstLow, firstHigh)] << endl;

	  //      cout << "\tat global distances from mean: " << prox[ Set::FullSet(firstLow, firstHigh)] << endl;

	  // check whether between partitions is even necessary
	  if ( sparsify.isZero((prox[secondLow] - prox[firstHigh])/2) )
	    {
	      // find between-set links and within set links (merge)
	      cout << "next block is close: " << prox[secondLow]-prox[firstHigh] << endl;
	      cout << "AxA" << endl;

	      linkSimilar( cat[ trueFromPoints(catFromBools, 1, firstLow, firstHigh) ] , cat[ trueFromPoints(catFromBools, 0, firstLow, firstHigh) ] );

	      cout << "AxB" << endl;
	      linkSimilar( cat[ trueFromPoints(catFromBools, 1, firstLow, firstHigh) ] , cat[ trueFromPoints(catFromBools, 0, secondLow, secondHigh) ] );

	      cout << "BxA" << endl;
	      linkSimilar( cat[ trueFromPoints(catFromBools, 1, secondLow, secondHigh) ] , cat[ trueFromPoints(catFromBools, 0, firstLow, firstHigh) ] );
	    }
	  else
	    {
	      cout << "next block is far: " << prox[secondLow]-prox[firstHigh] << endl;
	      // find within-set links for each
	      linkSimilar( cat[ trueFromPoints(catFromBools, 1, firstLow, firstHigh) ] , cat[ trueFromPoints(catFromBools, 0, firstLow, firstHigh) ] );
	    }

	  if ( k == beg.size()-2 )
	    {
	      //	  cout << "doing last block with itself" << endl;
	      linkSimilar( cat[ trueFromPoints(catFromBools, 1, secondLow, secondHigh) ] , cat[ trueFromPoints(catFromBools, 0, secondLow, secondHigh) ] );
	    }
	}
    }

  //  cout << "Done block link similar" << endl;
}

Array<int> KDMethod::trueFromPoints(const Array<int> & catFromBools, int matchVal, int lowInd, int highInd)
{
  Array<int> result;
  for (int k=lowInd; k<=highInd; k++)
    {
      if (catFromBools[k] == matchVal) 
	result.add(k);
    }

  return result;
}

Array<int> KDMethod::blockBeginnings(const Array<double> & sortedProx)
{
  //  cout << "Block beg" << endl;

  int ind = 0;

  Array<int> result(1, 0);
  
  for (int k=0; k<sortedProx.size(); k++)
    {
//      cout << sortedProx[k] << " ";

      if ( sparsify.isNonzero( (sortedProx[k] - sortedProx[ind])/2 ) )
	{
	  ind = k;
	  
	  result.add(k);
	}

//        if ( k > 5 )
// 	 {
// 	   cout << sortedProx[ sortedProx.size() - 3 ] << " ";
// 	   cout << sortedProx[ sortedProx.size() - 2 ] << " ";
// 	   cout << sortedProx[ sortedProx.size() - 1 ] << endl;

// 	   cout << endl << result << endl;
// 	   exit(0);
// 	 }

    }

  //  cout << "Done block beg" << endl;
  return result;
}

void KDMethod::indent(int i) const
{
  for (int k=0; k<i; k++)
    {
      cout << '\t';
    }
}

void KDMethod::linkSimilar(const Array<int> & from, const Array<int> & to, bool forced)
{
  recCount++;
  indent(recCount);
  cout << recCount << " " << "link " << from.size() << "x" << to.size() << endl;

  if ( ! forced && from.size() * to.size() > 10000 )
    {
      blockLinkSimilar(from, to);
    }
  else
    {
//       cout << "Linking similar:" << endl;
//       cout << "\t" << from << endl;
//       cout << "x\t" << to << endl << endl;

      if ( forced )
	cout << "Forcing " << from.size() <<"x" << to.size() << " = " << from.size()*to.size() << endl;

      for( int k=0; k<from.size(); k++)
	{
	  for ( int j=0; j<to.size(); j++)
	    {
	      if ( sparsify.isZero( diffNormSquared(*orderedFeatureVectorRefs[ from[k] ], *orderedFeatureVectorRefs[ to[j] ]) ) )
		{
		  similar[ from[k] ].add( to[j] );
		  //		  cout << "\t\tEdge " << from[k] << " --> " << to[j] << endl;
		}
	    }
	}
    }
  recCount--;
  //  cout << "Done low n^2 link" << endl;
}

Vector KDMethod::meanVector(const Array<int> & s)
{
  // assumes the array is a valid set; no element may appear more than once
  // also assumes that the array is not empty

  Vector res = *orderedFeatureVectorRefs[ s[0] ];
  
  for (int k=1; k<s.size(); k++)
    {
      res += *orderedFeatureVectorRefs[ s[k] ];
    }

  res /= s.size();

  return res;
}

Array<Array<int> > KDMethod::getSimilar()
{
  int upper = orderedFeatureVectorRefs.size() - 1;
  linkSimilar(Set::FullSet(0, upper), Set::FullSet(0, upper));

  for (int k=0; k<similar.size(); k++)
    {
      if ( similar[k].size() == 0 )
	cerr << "Error: vector " << k << " has no neighbors (it should have itself at least" << endl;
    }

  return similar;
}

Array<Array<int> > KDMethod::slowSimilar()
{
  for (int k=0; k<orderedFeatureVectorRefs.size(); k++)
    {
      for (int j=0; j<orderedFeatureVectorRefs.size(); j++)
	{
	  if ( sparsify.isZero( diffNormSquared(*orderedFeatureVectorRefs[ k ], *orderedFeatureVectorRefs[ j ]) ) )
	    {
	      //	      cout << "points " << k << " and " << j << " are this far " << diffNormSquared(*orderedFeatureVectorRefs[ k ], *orderedFeatureVectorRefs[ j ]) << endl;
	      similar[k].add(j);
	    }
	}
    }

  return similar;
}

