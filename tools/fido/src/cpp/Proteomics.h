#ifndef _PROTEOMICS_H
#define _PROTEOMICS_H

#include "Array.h"
#include <fstream>
#include "BasicBigraph.h"

using namespace std;

class Proteomics
{
 public:
  static Array<string> trypsinDigest(const string & str);
  static bool observablePeptide(const string & str);
  static void setMinMaxLength(char * pivdoFname, double minScore);
  static int maxLength, minLength;
};

#endif
