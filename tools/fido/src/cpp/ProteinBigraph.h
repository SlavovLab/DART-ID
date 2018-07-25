#ifndef _ProteinBigraph_H
#define _ProteinBigraph_H

#include "StringTable.h"
#include "Array.h"
#include "Set.h"
#include "Matrix.h"
#include "Random.h"

using namespace std;

class ProteinBigraph
{
public:
  ProteinBigraph();

  friend istream & operator >>(istream & is, ProteinBigraph & sb);
  void outputLP() const;

  void printProteinWeights() const;
  void sumProteins();

  void EM(int K);
  void scoreProteins();

  void partitionSections();

  class ProteinBigraphFormatException: public exception {
   virtual const char* what() const throw()
   {
    return "protein bigraph: format exception";
   }
  };
protected:

  struct GraphNode
  {
    const string & name;
    Set & association;
    double & weight;
    int & section;

  GraphNode(const string & n, Set & a, double & w, int & s):
    name(n),
      association(a),
      weight(w),
      section(s)
    {
    }
  };

  struct GraphLayer
  {
    StringTable names;
    Array<Set> associations;
    Array<double> weights;
    Array<int> sections;
    GraphLayer*partner;

    GraphNode operator [] (int k)
    {
      return GraphNode(names[k], associations[k], weights[k], sections[k]);
    }
    int size()
    {
      // all should be the same
      return associations.size();
    }
  };

  void read(istream & is);
  void add(GraphLayer & gl, const string & item);
  void connect(const string & pepStr, const string & protStr);

  void disconnectPSM(int k);
  void removePoorPSMs();
  void reindex();

  void traceConnected(int index, GraphLayer & gl, int sectionNumber);

  double alpha, beta;

  double Threshold;

  GraphLayer proteinsToPSMs, PSMsToProteins;
};

#endif

