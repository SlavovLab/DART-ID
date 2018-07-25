#ifndef _FASTAREADER_H
#define _FASTAREADER_H

#include "Array.h"
#include <fstream>

class FastaReader
{
 public:
  FastaReader(char*fname);
  void read(char*fname);
  void trimAsterisk(string & seq) const;
  Array<string> names;
  Array<string> sequences;
};


#endif
