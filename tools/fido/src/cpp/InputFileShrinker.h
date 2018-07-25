#ifndef _InputFileShrinker_H
#define _InputFileShrinker_H


#include <cstdlib>
#include <iostream>
#include <ios>
#include <istream>
#include <ostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <list>
#include <iomanip>

using namespace std;

class FormatException :public exception {
   virtual const char* what() const throw()
   {
    return "array: format exception";
   }
  };

class SizeException : public exception {
   virtual const char* what() const throw()
   {
    return "dataset: size exception";
   }
  };


class InputFileShrinker
{
public:

  InputFileShrinker();

  void shrink(const char * inputfilename, const char * outputfilename, unsigned int const targetSize);
 
private: 
  void shrink(istream & is, ostream & os, unsigned int const targetSize);
  void read(istream & is);
  void reduce(unsigned int const targetSize);
  void write(ostream & os);
  
  
  map<int, double> 				m_chargeStates;
  map<int, double>::iterator 	m_iterChargeStates;
  
  // pepname, prob, chargestate, set of proteins
  multimap<string, pair< pair<double, int>, set<string> > > m_data;
  multimap<string, pair< pair<double, int>, set<string> > >::iterator  m_iterdata;
  
};

#endif

