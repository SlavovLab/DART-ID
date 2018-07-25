#ifndef _GramSchmidt_H
#define _GramSchmidt_H

#include "Matrix.h"

class GramSchmidt : public Matrix
{
public:
  GramSchmidt()
    { }

 GramSchmidt(int c) :
  Matrix(c)
    {
    }

  Vector projectOrthogonal(const Vector & rhs) const;
  void add(const Vector & rhs);

  class ZeroVectorException : public exception {
   virtual const char* what() const throw()
   {
    return "GramSchmidt: zero vector exception";
   }
  };
  void swapToLast(int a);

protected:
  void swapWithNext(int k);

  Vector projectSingle(const Vector & basis, const Vector & v) const;

  Matrix original;
  static Numerical zeroVectorChecker;
};

#endif

