#ifndef _Set_H
#define _Set_H

#include "Array.h"

class Set : public Array<int>
{
 public:
  Set() {}
  Set(Array<int> unsorted);

  void add(int el);
  const Set & operator &=(const Set & rhs);
  const Set & operator |=(const Set & rhs);
  bool isEmpty() const
  {
    return size() == 0;
  }

  int find(int x) const;
  int findHelper(int low, int high, int x) const;

  Set reindexTo(const Set & rhs);
  Set reindexToFind(const Set & rhs);

  class UnsortedOrderException: public exception {
   virtual const char* what() const throw()
   {
    return "set: unsorted order exception";
   }};
  class InvalidBaseException :public exception {
   virtual const char* what() const throw()
   {
    return "set: invalid base exception";
   }
  };


  Array<int> operator [](const Set & rhs) const
  {
    return Array<int>::operator [](rhs);
  }

  const int & operator [](int k) const
  {
    return Array<int>::operator [](k);
  }

  // for hashing sets
  static unsigned int sumSetElements(const Set & s)
  {
    unsigned int sum = 0;
    for (int k=0; k<s.size(); k++)
      sum += s[k];

    return sum;
  }

  using Array<int>::size;
  using Array<int>::Iterator;
  using Array<int>::begin;
  using Array<int>::end;

  Set without( const Set & rhs ) const;

  static Set FullSet(int low, int high);
  static Set SingletonSet(int value)
  {
    return FullSet(value, value);
  }
  
  int randomElement() const
  {
    return (*this)[ rand() % size() ];
  }

  bool contains(int ind) const;

 private:

  bool verify() const;
};

Set operator &(Set lhs, const Set & rhs);
Set operator |(Set lhs, const Set & rhs);

Array<int> inverseMap(const Array<int> & forwardMap);

//ostream & operator <<(ostream & os, const Set & rhs);
//istream & operator >>(istream & is, Set & rhs);

#endif

