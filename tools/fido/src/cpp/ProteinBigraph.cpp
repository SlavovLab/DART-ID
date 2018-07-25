#include "ProteinBigraph.h"

ProteinBigraph::ProteinBigraph()
{
  Threshold = .2;
  //  Threshold = 0.0;
  proteinsToPSMs.partner = & PSMsToProteins;
  PSMsToProteins.partner = & proteinsToPSMs;
}

istream & operator >>(istream & is, ProteinBigraph & sb)
{
  sb.read(is);
  return is;
}

void ProteinBigraph::read(istream & is)
{
  char instr;
  string pepName, protName;
  double value;

  int pepIndex = -1;

  while ( is >> instr )
    {
      if ( instr == 'e' )
	{
	  is >> pepName;

	  if ( PSMsToProteins.names.lookup(pepName) == -1 )
	    add(PSMsToProteins, pepName);
	  
	  pepIndex = PSMsToProteins.names.lookup(pepName);
	}
      else if ( instr == 'r' )
	{
	  is >> protName;

	  if ( proteinsToPSMs.names.lookup(protName) == -1 )
	    add(proteinsToPSMs, protName);

	  // this is a hack to only get the nodes connected the first
	  // time the peptide is seen; (it's initial weight is -1.0
	  // and will always be changed to >= 0.0 when it is first
	  // finished being seen. Might want to make this better
	  // later.

	  if ( PSMsToProteins.weights[ pepIndex ] == -1.0 )
	    connect(pepName, protName);
	}
      else if ( instr == 'p' )
	{
	  is >> value;

	  PSMsToProteins.weights[ pepIndex ] = max(PSMsToProteins.weights[pepIndex], value);
	}
      else
	{
	  throw ProteinBigraphFormatException();
	}
    }

  /***
  cout << "# PSMs\t" << PSMsToProteins.size() << endl;
  cout << "\t" << PSMsToProteins.weights << endl;
  cout << "PSM  assoc\t" << PSMsToProteins.associations << endl;
  cout << "Prot assoc\t" << proteinsToPSMs.associations << endl;
  ***/

  cout << "Read, removing poor" << endl;
  removePoorPSMs();

  cout << "Reindexing" << endl;
  reindex();

  /***
  cout << "# PSMs\t" << PSMsToProteins.size() << endl;
  cout << "\t" << PSMsToProteins.weights << endl;
  cout << "PSM  assoc\t" << PSMsToProteins.associations << endl;
  cout << "Prot assoc\t" << proteinsToPSMs.associations << endl;
  ***/

  cout << "Constructed" << endl;
}

void ProteinBigraph::reindex()
{
  GraphLayer newPSMs, newProteins;

  int k;

  // note: later check to see if you can do suchThat in Array using a
  // pointer to a member function
  Set connectedPSMs;
  for (k=0; k<PSMsToProteins.size(); k++)
    {
      if ( ! PSMsToProteins.associations[k].isEmpty() )
	connectedPSMs.add(k);
    }

  Set connectedProteins;
  for (k=0; k<proteinsToPSMs.size(); k++)
    {
      if ( ! proteinsToPSMs.associations[k].isEmpty() )
	connectedProteins.add(k);
    }

  cout << "Connected found" << endl;

  cout << "\tnames" << endl;
  newPSMs.names = StringTable::AddElements( PSMsToProteins.names.getItemsByNumber()[ connectedPSMs ] );
  cout << "\tassociations" << endl;
  newPSMs.associations = PSMsToProteins.associations[ connectedPSMs ];
  cout << "\tweights" << endl;
  newPSMs.weights = PSMsToProteins.weights[ connectedPSMs ];
  cout << "\tsections" << endl;
  newPSMs.sections = PSMsToProteins.sections[ connectedPSMs ];

  cout << "\tReindexing association sets" << endl;

  for (k=0; k<newPSMs.associations.size(); k++)
    {
      newPSMs.associations[k] = newPSMs.associations[k].reindexTo( connectedProteins );
    }

  cout << "Done PSMs" << endl;

  cout << "\tnames" << endl;
  newProteins.names = StringTable::AddElements( proteinsToPSMs.names.getItemsByNumber()[ connectedProteins ] );
  cout << "\tnames" << endl;
  newProteins.associations = proteinsToPSMs.associations[ connectedProteins ];
  cout << "\tnames" << endl;
  newProteins.weights = proteinsToPSMs.weights[ connectedProteins ];
  cout << "\tnames" << endl;
  newProteins.sections = proteinsToPSMs.sections[ connectedProteins ];
  cout << "\tnames" << endl;

  for (k=0; k<newProteins.associations.size(); k++)
    {
      newProteins.associations[k] = newProteins.associations[k].reindexTo( connectedPSMs );
    }

  cout << "Done proteins" << endl;
  
  proteinsToPSMs = newProteins;
  PSMsToProteins = newPSMs;

  proteinsToPSMs.partner = & PSMsToProteins;
  PSMsToProteins.partner = & proteinsToPSMs;
}

void ProteinBigraph::removePoorPSMs()
{
  int k;
  for (k=0; k<PSMsToProteins.size(); k++)
    {
      if ( PSMsToProteins.weights[k] < ProteinBigraph::Threshold )
	disconnectPSM(k);
    }
}

void ProteinBigraph::disconnectPSM(int k)
{
  Set & as = PSMsToProteins.associations[k];
  for (Set::Iterator iter = as.begin(); iter != as.end(); iter++)
    {
      Set & setRef = proteinsToPSMs.associations[ *iter ];
      setRef = setRef.without( Set::SingletonSet( k ) );
    }
  as = Set();
}

void ProteinBigraph::add(GraphLayer & gl, const string & item)
{
  if ( gl.names.lookup(item) == -1 )
    {
      // if the string is not already known, then add a new node for it
      gl.names.add(item);
      gl.associations.add( Set() );
      gl.weights.add( -1.0 );
      gl.sections.add(-1);
    }
}

void ProteinBigraph::connect(const string & pepName, const string & protName)
{
  int pepIndex = PSMsToProteins.names.lookup(pepName);
  int protIndex = proteinsToPSMs.names.lookup(protName);

  // performance note: currently O(n^2) worstcase. Later use a bitset
  // and then after the graph is read, pack it into a set. 

  //  cerr << "Joining " << protName << " and " << pepName << endl;

  PSMsToProteins.associations[ pepIndex ] |= Set::SingletonSet(protIndex);
  proteinsToPSMs.associations[ protIndex ] |= Set::SingletonSet(pepIndex);
}

/***
void ProteinBigraph::outputLP() const
{
  int n = proteinNames.numElements();
  cout << n << endl;

  cout << proteinNames.getItemsByNumber() << endl;

  Matrix A(n);
  Array<double> b;

  int k,j;

  // constraints for sum of proteins for a peptide <=1 
  for (k=0; k<PSMNodes.size(); k++)
    {
      Array<double> aRow(n, 0.0);

      const Set & set = PSMNodes[k].connectionSet;
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

  Vector one = Array<double>(PSMNodes.size(), 1.0);

  // make the objective vector
  Vector f(n);
  Array<double> fRow(n, 0.0);
  for (k=0; k<PSMNodes.size(); k++)
    {
      double pepWeight = PSMNodes[k].weight;

      const Set & s = PSMNodes[k].connectionSet;
      for (Set::Iterator si = s.begin(); si != s.end(); si++)
	{
	  //	  cerr << proteinNames.getItemsByNumber()[ *si ] << " gets " << pepWeight << endl;
	  fRow[ *si ] += pepWeight;
	}

      //      f += PSMNodes[k].weight * one[ PSMNodes[k].connectionSet ];


      //      cout << "Psst... PSM " << PSMNames.getItemsByNumber()[k] << " likes proteins: " << proteinNames.getItemsByNumber()[ PSMNodes[k].connectionSet ] << endl;
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
***/

 /***
void ProteinBigraph::sumProteins()
{
  for (int k=0; k<proteinNodes.size(); k++)
    {
      Node & pn = proteinNodes[k];

      //      cout << "Tabulating " << proteinNames.getItemsByNumber()[k] << endl;

      double sum = 0;
      for (Set::Iterator iter = pn.connectionSet.begin(); iter != pn.connectionSet.end(); iter++)
	{
	  Node & psmn = PSMNodes[ *iter ];

	  //	  cout << "\tpep " << psmn.weight << " -> " << (1-psmn.weight)+psmn.weight.connectionSet.size()-1/

	  sum += psmn.weight;

	  //	  sum -= 1/(psmn.weight*log(psmn.weight));
	}

      pn.weight = sum;
    }
}

void ProteinBigraph::scoreProteins()
{
  int k,j;

  Array<double> prodProtProbsForPSM(PSMNodes.size(), 0.0);
  Array<double> sumProtProbsForPSM(PSMNodes.size(), 0.0);
  for (j=0; j<prodProtProbsForPSM.size(); j++)
    {
      const Set & conSet = PSMNodes[j].connectionSet;

      double prod = 1.0;
      double tot = 0.0;
      for (Set::Iterator iter = conSet.begin(); iter != conSet.end(); iter++)
	{
	  prod *= (1-proteinNodes[ *iter ].weight) + proteinNodes[ *iter ].weight * (1-alpha);
	  tot += proteinNodes[ *iter ].weight;
	}

      sumProtProbsForPSM[j] = tot;
      prodProtProbsForPSM[j] = prod * (1-beta);
    }

  for (k=0; k<proteinNodes.size(); k++)
    {
      Node & pn = proteinNodes[k];

      //      cout << "Tabulating " << proteinNames.getItemsByNumber()[k] << endl;

      double prod = 1.0;

      for (Set::Iterator iter = pn.connectionSet.begin(); iter != pn.connectionSet.end(); iter++)
	{
	  Node & psmn = PSMNodes[ *iter ];

	  //	  cout << "\tpep " << psmn.weight << " -> " << (1-psmn.weight)+psmn.weight.connectionSet.size()-1/
	  
	  double denom = (1 - pn.weight) + pn.weight*(1-alpha);
	  if ( Vector::comparator.isNonzero( denom ) )
	    {
	      double term = pow( (1-psmn.weight) + psmn.weight * (1 - prodProtProbsForPSM[ *iter ] / denom) , sumProtProbsForPSM[ *iter ] * (alpha) );

	      if ( isnan(term) )
		{
		  cerr << "Nan found. " << pn.weight << " " << prodProtProbsForPSM[ *iter ] << " " << sumProtProbsForPSM[ *iter ] << endl;
		}
	      prod *= term;

	    }
	}

      pn.weight = 1 - prod;
    }
}
***/

/***
void ProteinBigraph::EM(int K)
{
  int rep;
  for (rep = 0; rep<K; rep++)
    {
      scoreProteins();
 
      Array<double> rArray( proteinNodes.size() );

      //  cout << r << endl << r.unpack() << endl;

      for (int k=0; k<proteinNodes.size(); k++)
	{
	  rArray[k] = proteinNodes[k].weight;
	}
      Vector r = rArray;

      //      alpha = estimateAlpha(r);
      //      beta = estimateBeta(r);
      
      //      cerr << "REP " << rep << endl << endl;
      //      printProteinWeights();
      //      cerr << endl;

      //      scorePeptides();
      if ( rep % 100 == 0 )
	{
	  //	  cout << "REP " << rep << endl << endl;
	  //	  printProteinWeights();
	}
    }
  //  cout << "REP " << rep << endl << endl;
  //  printProteinWeights();
}
***/

/***
void ProteinBigraph::printProteinWeights() const
{
  vector<pair<double, string> > sortie(proteinNodes.size());
  const Array<string> & prots = proteinNames.getItemsByNumber();

  int k;
  for (k=0; k<proteinNodes.size(); k++)
    {
      //     cerr << proteinNodes[k].weight << "\t" << prots[k] << endl;
      sortie[k] = pair<double, string>(proteinNodes[k].weight, prots[k]);
    }

  sort( sortie.begin(), sortie.end() , greater<pair<double, string> >() );
  //  cerr << "Sorting" << endl;

  for (k=0; k< min(300, proteinNodes.size() ); k++)
    //    for (k=0; k< proteinNodes.size(); k++)
    {
      // trim out 0s for speed?
      //      if ( sortie[k].first > 0.0 )
      cout << sortie[k].first << " { " << sortie[k].second << " } " << endl;
    }
}
***/

void ProteinBigraph::traceConnected(int index, GraphLayer & gl, int sectionNumber)
{
  if ( gl.sections[index] == sectionNumber )
    return;

  gl.sections[index] = sectionNumber;

  const Set & as = gl.associations[index];
  for (Set::Iterator iter = as.begin(); iter != as.end(); iter++)
    {
      traceConnected( *iter, *gl.partner, sectionNumber );
    }
}

void ProteinBigraph::partitionSections()
{
  int section = 0;
  int k;
  cout << "Tracing proteins" << endl;
  for (k=0; k<proteinsToPSMs.size(); k++)
    {
      if ( proteinsToPSMs[k].section == -1 )
	{
	  traceConnected( k, proteinsToPSMs, section );
	  section++;
	}
    }

  cout << "Counting sizes of sections" << endl;
  Array<int> proteinSectionSizes(section, 0), PSMSectionSizes(section, 0);
  for (k=0; k<proteinsToPSMs.size(); k++)
    {
      proteinSectionSizes[ proteinsToPSMs[k].section ]++;
    }

  for (k=0; k<PSMsToProteins.size(); k++)
    {
      PSMSectionSizes[ PSMsToProteins[k].section ]++;
    }

  cout << "# sections = " << section << endl;
  cout << "\tproteins :" << proteinSectionSizes << endl;
  cout << "\tPSMs     :" << PSMSectionSizes << endl;
}

