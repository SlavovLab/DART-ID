#include "FastaReader.h"

FastaReader::FastaReader(char*fname)
{
  read(fname);
}

void FastaReader::trimAsterisk(string & seq) const
{
  if ( seq.length() > 0 && seq[ seq.length()-1 ] == '*' )
    seq.resize(seq.length()-1);
}

void FastaReader::read(char*fname)
{
  names = Array<string>();
  sequences = Array<string>();

  ifstream fin(fname);
  string buff;

  string name = "", seq = "";

  while (getline(fin, buff))
    {
      if ( buff == "" )
	continue;

      istringstream ist(buff);
      char c;
      ist >> c;

      if ( c == '>' )
	{
	  // do not add the first one
	  if ( name != "" )
	    {
	      trimAsterisk(seq);
	      sequences.add(seq);
	    }

	  ist >> name;
	  names.add(name);
	  seq = "";
	  continue;
	}
      else
	{
	  seq += buff;
	}
    }
  if ( name != "" )
    {
      trimAsterisk(seq);
      sequences.add(seq);
    }
}

