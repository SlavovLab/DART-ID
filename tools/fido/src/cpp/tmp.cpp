#include "Cache.h"
#include <iostream>
#include <vector>
#include "Model.h"

using namespace std;

ostream & operator <<(ostream & os, const vector<double> & rhs)
{
  os << "{ ";
  for (int k=0; k<rhs.size(); k++)
    {
      os << rhs[k] << " \t";
    }
  os << "}";
  return os;
}

class Terminal
{
protected:
  int size;
public:
  Terminal(int s)
  {
    size = s;
  }
  vector<double> sequence(double scale) const
  {
    cout << "Calling sequence" << endl;
    vector<double> result;

    for (int k=0; k<=size; k++)
      {
	result.push_back(k*scale);
      }

    return result;
  }
  double seqSum(double scale) const
  {
    cout << "Calling seqSum" << endl;
    vector<double> s = sequence(scale);

    double tot = 0.0;

    for (int k=0; k<s.size(); k++)
      {
	tot += s[k];
      }

    return tot;
  }
};

int main()
{
  Terminal term(4);
  //  MemberFunction<Terminal, vector<double>, double> sdb( &Terminal::sequence, &term);
  
  Model a(.1, .1, .2);
  Model b(.1, .2, .1);

  if ( a ==  b )
    cout << "Equal" << endl;
  else
    cout << "Nequal" << endl;

  return 0;
}

