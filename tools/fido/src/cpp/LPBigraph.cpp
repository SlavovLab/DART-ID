#include "LPBigraph.h"

void LPBigraph::outputLP() const
{
  int n = proteinsToPSMs.size();
  cout << n << endl;

  cout << proteinsToPSMs.names.getItemsByNumber() << endl;

  Matrix A(n);
  Array<double> b;

  int k,j;

  // constraints for sum of proteins for a peptide <=1 
  for (k=0; k<PSMsToProteins.size(); k++)
    {
      Array<double> aRow(n, 0.0);

      const Set & set = PSMsToProteins[k].association;
      for (Set::Iterator setIter = set.begin(); setIter != set.end(); setIter++)
	{
	  aRow[ *setIter ] = 1.0;
	}

      if ( norm(aRow) > 0.0 )
	{
	  A.add( -Vector(aRow) );
	  b.add( -1.0 );
	}
    }

  // constraints for proteins <= 1
  for (j=0; j<n; j++)
    {
      Vector protConstraint(n);
      protConstraint.addElement(j, 1.0);

      A.add( -protConstraint );
      b.add( -1.0 );
    } 

  // constraints for proteins >= 0
  for (j=0; j<n; j++)
    {
      Vector protConstraint(n);
      protConstraint.addElement(j, 1.0);

      A.add( protConstraint );
      b.add( 0.0 );
    } 

  // constraint for spending few proteins
  A.add( Vector( Array<double>(n, -1.0) ) );
  // initial gamma
  b.add( -0.0 );

  Vector one = Array<double>(PSMsToProteins.size(), 1.0);

  // make the objective vector
  Vector f(n);
  Array<double> fRow(n, 0.0);
  for (k=0; k<PSMsToProteins.size(); k++)
    {
      double pepWeight = PSMsToProteins[k].weight;

      const Set & s = PSMsToProteins[k].association;
      for (Set::Iterator si = s.begin(); si != s.end(); si++)
	{
	  //	  cerr << proteinNames.getItemsByNumber()[ *si ] << " gets " << pepWeight << endl;
	  fRow[ *si ] += pepWeight;
	}

      //      f += PSMNodes[k].weight * one[ PSMNodes[k].association ];


      //      cout << "Psst... PSM " << PSMNames.getItemsByNumber()[k] << " likes proteins: " << proteinNames.getItemsByNumber()[ PSMNodes[k].association ] << endl;
      //      cout << "\tand has weight " << PSMNodes[k].weight << endl << endl;
      //      cout << "\tand index " << k << endl;
      //      getchar();
    }


  f = fRow;

  Vector x(n);
  cout << x << endl;
  cout << f << endl;
  cout << Vector(b) << endl;
  cout << A << endl;
}

void LPBigraph::sumProteins()
{
  for (int k=0; k<proteinsToPSMs.size(); k++)
    {
      Node & pn = proteinsToPSMs[k];

      //      cout << "Tabulating " << proteinNames.getItemsByNumber()[k] << endl;

      double sum = 0;
      for (Set::Iterator iter = pn.association.begin(); iter != pn.association.end(); iter++)
	{
	  Node & psmn = PSMsToProteins[ *iter ];

	  //	  cout << "\tpep " << psmn.weight << " -> " << (1-psmn.weight)+psmn.weight.association.size()-1/

	  sum += psmn.weight;

	  //	  sum -= 1/(psmn.weight*log(psmn.weight));
	}

      pn.weight = sum;
    }
}
