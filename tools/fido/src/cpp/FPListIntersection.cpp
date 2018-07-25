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

int intersectionSize(Array<string> A, Array<string> B)
{
  A.sort();
  B.sort();

  //  cout << "\tafter sorting" << endl;
  //  cout << "\t\t" << A << endl << "\t\t" << B << endl;

  int c=0;

  for (int kA=0, kB=0; kA<A.size() && kB<B.size(); )
    {
      //      cout << "\t\t\tcurrent is: " << A[kA] << " " << B[kB] << endl;
      if ( A[kA] == B[kB] )
	{
	  //	  cout << "\t\toverlap: " << A[kA] << " " << B[kB] << endl;
	  c++;

	  kA++;
	  kB++;
	}
      else
	{
	  if ( A[kA] > B[kB] )
	    kA++;
	  else
	    kB++;
	}
    }

  //  cout << "Overlap size was: " << c << endl;

  return c;
}

int main(int argc, char**argv)
{
  if ( argc == 3 )
    {
      ifstream finA(argv[1]), finB(argv[2]);

      Array<int> fpA, fpB;
      Array<Array<string> > tpA, tpB;

      finA >> fpA >> tpA;
      finB >> fpB >> tpB;

      // this code assumes that one of the arrays has at least one element
      int x = 0 , y = 0;
      int kA, kB;
      int f;
      for (kA=0, kB=0; ; )
	{
	  //	  cout << "On " << fpA[kA] << "," << tpA[kA].size() << "        " << fpB[kB] << "," << tpB[kB].size() << endl;

	  f = max( fpA[kA], fpB[kB] );

	  //	  if ( f > x && ( kA != 0 || kB != 0 ) )
	    {
	      // ready to update
	      cout << x << " " << y << endl;
	    }

	  x=f;
	  y=intersectionSize(tpA[kA], tpB[kB]);

	  // stopping criteria
	  if ( kA == fpA.size()-1 && kB == fpB.size()-1 )
	    break;

	  if ( kA == fpA.size()-1 || kB == fpB.size()-1 )
	    {
	      //	      cout << "Edge advance" << endl;

	      if ( kA == fpA.size()-1 )
		kB++;
	      else
		kA++;
	    }
	  else
	    { 
	      // advance the correct one
	      if ( fpA[kA+1] < fpB[kB+1] )
		{
		  //		  cout << "Advance A" << endl;
		  kA++;
		}
	      else if ( fpA[kA+1] > fpB[kB+1] )
		{
		  //		  cout << "Advance B" << endl;
		  kB++;
		}
	      else
		{
		  //		  cout << "Advance both" << endl;
		  kA++;
		  kB++;
		}
	    }
	}

      //      if ( f > x && ( kA != 0 || kB != 0 ) )
	{
	  // tail update
	  cout << x << " " << y << endl;
	}
    }
  else
    {
      cerr << "usage: <tp-file 1> <tp-file 2>" << endl << endl;
    }
  return 0;
}


