#include "Numerical.h"

bool Numerical::isPos(double d) const
{
  return d > epsilon;
}

bool Numerical::isNonpos(double d) const
{
  return d <= epsilon;
}

bool Numerical::isNeg(double d) const
{
  return d < -epsilon;
}

bool Numerical::isNonneg(double d) const
{
  return d >= -epsilon;
}

bool Numerical::isZero(double d) const
{
  return fabs(d) <= epsilon;
}

bool Numerical::isNonzero(double d) const
{
  return fabs(d) > epsilon;
}

bool Numerical::isEqual(double a, double b) const
{
  return isZero( b - a );
}

bool Numerical::isInequal(double a, double b) const
{
  return isNonzero( b - a );
}

bool Numerical::isDifferentSign(double a, double b) const
{
  return isPos(a) && isNeg(b) || isNeg(a) && isPos(b);
}
