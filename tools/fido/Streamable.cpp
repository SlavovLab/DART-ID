#include "Streamable.h"

ostream & operator << (ostream & os, const Streamable & rhs)
{
  return rhs.putToStream(os);
}

istream & operator >> (istream & is, Streamable & rhs)
{
  return rhs.getFromStream(is);
}

