#include "Array.h"


Set::Set(Array<int> unsorted)
{
  unsorted.sortA();
  
  for (int k=0; k<unsorted.size(); k++)
    {
      add(unsorted[k]);
    }
}

void Set::add(int el)
{
  #ifdef SAFE_ARRAYS
  if ( ! isEmpty() && el <= back() )
    {
      throw UnsortedOrderException();
    }
  #endif

  Array<int>::add(el);
}

const Set & Set::operator &=(const Set & rhs)
{
  Set result = (*this) & rhs;
  return *this = result;
}

int Set::find(int x) const
{
  return Array<int>::sortedFind(x);
}

int Set::findHelper(int low, int high, int x) const
{
  return Array<int>::sortedFindHelper(low, high, x);
}

Set Set::reindexTo(const Set & base)
{
  //  C == A.reindexTo(B) iff 
  //  B[C] == A

  Set C;

  if ( base.isEmpty() )
    throw InvalidBaseException();

  for (Set::Iterator iter = begin(), baseIter = base.begin(); iter != end() && baseIter != base.end(); )
    {
      if ( *iter > *baseIter )
	{
	  baseIter++;
	}
      else if ( *iter == *baseIter )
	{
	  C.add( baseIter.getLocation() );
	  iter++;
	  baseIter++;
	}
      else
	{
	  // skips an element of iter
	  throw InvalidBaseException();
	}
    }

  return C;
}

Set Set::reindexToFind(const Set & base)
{
  //  C == A.reindexTo(B) iff 
  //  B[C] == A

  Set C;

  if ( base.isEmpty() )
    {
      cerr << "Empty base in reindex-- set size is " << size() << endl;
      throw InvalidBaseException();
    }

  int lastFind = 0;

  for (Set::Iterator iter = begin(); iter != end(); iter++)
    {
      int loc = base.findHelper(lastFind, base.size()-1, *iter);

      //      cout << "\t\t\tOld school " << (base == *iter) << endl;
      
      if ( loc == -1 )
	{
	  // skips an element of iter
	  cerr << "skipped in reindex-- set size is " << size() << endl;
	  throw InvalidBaseException();
	}

      //      cout << "\t\t\tfound " << loc << endl;

      lastFind = loc;

      C.add(loc);
    }

  return C;
}

const Set & Set::operator |=(const Set & rhs)
{
  Set result = (*this) | rhs;
  return *this = result;
}

Set operator &(const Set lhs, const Set & rhs)
{
  Set result;
  for (Set::Iterator iterLhs = lhs.begin(), iterRhs = rhs.begin(); iterLhs != lhs.end() && iterRhs != rhs.end(); )
    {
      if ( *iterLhs == *iterRhs )
	{
	  result.add(*iterLhs);
	  iterLhs++;
	  iterRhs++;
	}
      else if ( *iterLhs < *iterRhs )
	{
	  iterLhs++;
	}
      else
	iterRhs++;
    }

  return result;
}

Set operator |(const Set lhs, const Set & rhs)
{
  Set result;
  Set::Iterator iterLhs = lhs.begin(), iterRhs = rhs.begin();
  for ( ; iterLhs != lhs.end() && iterRhs != rhs.end(); )
    {
      if ( *iterLhs == *iterRhs )
	{
	  result.add(*iterLhs);
	  iterLhs++;
	  iterRhs++;
	}
      else if ( *iterLhs < *iterRhs )
	{
	  result.add(*iterLhs);
	  iterLhs++;
	}
      else
	{
	  result.add(*iterRhs);
	  iterRhs++;
	}
    }
  
  for ( ; iterLhs != lhs.end(); iterLhs++ )
    {
      result.add( *iterLhs );
    }
  for ( ; iterRhs != rhs.end(); iterRhs++ )
    {
      result.add( *iterRhs );
    }

  return result;
}

bool Set::verify() const
{
  for (int k=0; k<size()-1; k++)
    {
      if ( (*this)[k] >= (*this)[k+1] )
	throw UnsortedOrderException();
    }

  return true;
}

/***
istream & Set::getFromStream(istream & is)
{
  char delim;
  is >> delim;

  clear();

  // the character for an empty array
  if ( delim == '0')
    return is;

  if ( delim != '{' )
    {
      throw FormatException();
    }
  int element;
  
  // now that a { character has been read
  // it is guaranteed that at least one
  // element will occur

  do
    {
      is >> element;
      add(element);
      is >> delim;
      
      if ( delim != ',' && delim != '}')
	throw FormatException();

    } while ( delim != '}' );
  return is;
}

ostream & operator <<(ostream & os, const Set & rhs)
{
  return rhs.putToStream(os);
}

istream & operator >>(istream & is, Set & rhs)
{
  return rhs.getFromStream(is);
}
***/

Set Set::without( const Set & rhs ) const
{
  Set result;

  Set::Iterator iterThis, iterRhs;
  for (iterThis = begin(), iterRhs = rhs.begin(); iterThis != end() && iterRhs != rhs.end(); )
    {
      //      cout << "Rep : " << *iterThis << " " << *iterRhs << endl;
      if ( *iterThis < *iterRhs )
	{
	  result.add( *iterThis );
	  iterThis++;
	}
      else if ( *iterThis > *iterRhs )
	{
	  // do not add elements from rhs
	  iterRhs++;
	}
      else
	{
	  // they are equal, so do not add
	  iterThis++;
	  iterRhs++;
	}
    }

  for ( ; iterThis != end(); iterThis++)
    {
      //      cout << "Post rep : " << *iterThis << endl;
      result.add( *iterThis );
    }

  return result;
}

Set Set::FullSet(int low, int high)
{
  Set result;
  for (int k=low; k<=high; k++)
    {
      result.add(k);
    }

  return result;
}

Array<int> inverseMap(const Array<int> & forwardMap)
{
  Array<int> result(forwardMap.size(), -1);
  
  int k;
  for (k=0; k<result.size(); k++)
    {
      int ind = forwardMap[k];
      if ( result[ind] == -1 )
	result[ ind ] = k;
      else
	{
	  cerr << "Warning: inverseMap called without an injective forward map" << endl << endl;
	  exit(1);
	}
    }

  // check for -1s
  for (k=0; k<result.size(); k++)
    {
      if ( result[k] == -1 )
	{
	  cerr << "Warning: inverseMap called without an onto forward map" << endl << endl;
	  exit(1);
	}
    }
  
  return result;
}

bool Set::contains(int ind) const
{
  return find(ind) != -1;
}


