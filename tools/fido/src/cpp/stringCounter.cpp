#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include "Vector.h"
#include "StringTable.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 2 )
    {
      Array<Array<string> > allLists;

      ifstream fin(argv[1]);

      string buff;
      while ( getline(fin, buff) )
	{
	  istringstream ist(buff);
	  Array<string> list;
	  ist >> list;

	  allLists.add(list);
	}
      
      cout << "Made list of lists. Making StringTable" << endl;

      StringTable st(100001);
      
      int k, j;

      Array<int> counts;
      Array<string> names;
      for (k=0; k<allLists.size(); k++)
	{
	  const Array<string> & line = allLists[k];

	  for (j=0; j<line.size(); j++)
	    {
	      if ( st.add(line[j]) )
		{
		  counts.add(1);
		  names.add( line[j] );
		}
	      else
		{
		  counts[ st.lookup(line[j]) ]++;
		}
	    }
	}

      vector<pair<int, string> > sortie(counts.size());

      for (k=0; k<counts.size(); k++)
	{
	  sortie[k] = pair<int, string>(counts[k], names[k]);
	}

      for (k=0; k<counts.size(); k++)
	{
	  cout << sortie[k].first << "\t" << sortie[k].second << endl;
	}

      sort( sortie.begin(), sortie.end() );

      for (k=0; k<counts.size(); k++)
	{
	  cout << sortie[k].first << "\t" << sortie[k].second << endl;
	}
    }
  else
    {
      cerr << "usage: StringCounter <multi lists file>" << endl;
    }
  return 0;
}


