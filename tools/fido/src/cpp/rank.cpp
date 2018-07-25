#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include "Vector.h"

using namespace std;

int findRankBound(const Array<double> & a, int k1)
{
  int k2;
  for (k2 = k1; k2<a.size()-1; k2++)
    {
      if ( a[k2+1] != a[k1] )
	break;
    }

  // there are two cases arriving here:
  // 1. break statement (trivially OK)
  // 2. k2 = a.size()-1 and all between k1, k2 are equal

  return k2;
}

Array<double> getRanking(const Array<double> & a)
{
  Array<double> result(a.size(), -1);

  double window = 0;
  for (int k1 = 0; k1 < a.size(); )
    {
      int k2 = findRankBound(a, k1);

      double windowSize = k2-k1 + 1;
      double averageRank = ( ( window + windowSize )*( window + windowSize + 1 ) - window * (window + 1) ) / ( 2*windowSize );
      //      cout << "Giving window " << k1 << ", " << k2 << "\t" << averageRank << endl;

      for (int k=k1; k<=k2; k++)
	{
	  result[k] = averageRank / a.size();
	}
      
      window += k2 - k1 + 1;
      k1 = k2+1;
    }

  return result;
}

int main(int argc, char**argv)
{
  if ( argc == 2 )
    {
      string ascDes = argv[1];

      Array<double> a, backup;

      double el;
      while ( cin >> el )
	{
	  a.add(el);
	}

      if ( a.size() == 0 )
	{
	  cerr << "Error: cannot rank an array of size 0" << endl << endl;
	  exit(1);
	}

      backup = a;
      Array<int> order;
      
      if ( ascDes == "A" || ascDes == "a" )
	order = a.sortA();
      else if ( ascDes == "D" || ascDes == "d" )
	order = a.sort();
      else
	{
	  cerr << "Error: argument must be '+' or '-'" << endl << endl;
	  exit(1);
	}

      Array<double> ranking = getRanking(a);

      // now reorder to unsorted position (effectively a transpose)
      Array<int> unorder(order.size());

      int k;
      for (k=0; k<order.size(); k++)
	{
	  unorder[ order[k] ] = k;
	}

      for (k=0; k<unorder.size(); k++)
	{
	  cout << ranking[ unorder[k] ] << " " << backup[k] << endl;
	}
    }
  else
    {
      cerr << "usage: (for ascending order -- low values get low ranks ) Rank A < real_column_file" << endl << endl;
      cerr << "usage: (for descending order -- high values get low ranks ) Rank A < real_column_file" << endl << endl;
    }
  return 0;
}


