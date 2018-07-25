#include "InputFileShrinker.h"


InputFileShrinker::InputFileShrinker()
{
  
}



void InputFileShrinker::shrink(const char * inputfilename, const char * outputfilename, unsigned int const targetSize)
{
	// create two streams for that
	ifstream fin(inputfilename);
	ofstream fout(outputfilename);
	shrink(fin, fout, targetSize);
}

void InputFileShrinker::shrink(istream & is, ostream & os, unsigned int const targetSize)
{
	read(is);
	reduce(targetSize);
	write(os);
}

void InputFileShrinker::reduce(unsigned int const targetSize)
{
	if(targetSize >= m_data.size())
		return;
	
	
	if(targetSize < 51)
	{
		throw(SizeException());
	}
	
	multimap<string, pair< pair<double, int>, set<string> > > m_reduceddata;
	
	while(m_reduceddata.size() < targetSize)
	{
		// randomly remove an item
		m_iterdata = m_data.begin();
		std::advance( m_iterdata, rand()%m_data.size() );
		m_reduceddata.insert(*m_iterdata);
		m_data.erase(m_iterdata);
		if(m_reduceddata.size() % 263 == 0){
			cerr << "collected " << m_reduceddata.size() << " peptides for reduced parameter calculation" << endl;
		}
	}
	m_data = m_reduceddata;
	
}

void InputFileShrinker::write(ostream & os)
{
	os << "# this intermediate file was written by FIDO for accelarated target parameter estimation" << endl;
	
	// write the charge states
	
	os << setiosflags(ios::fixed) << setprecision(3);
	for(m_iterChargeStates = m_chargeStates.begin();m_iterChargeStates != m_chargeStates.end(); ++ m_iterChargeStates)
	{
		os << "d " << m_iterChargeStates->first << " " << m_iterChargeStates->second << endl;
	}
	
	os << setiosflags(ios::fixed) << setprecision(4);
	
	for(m_iterdata = m_data.begin(); m_iterdata != m_data.end(); ++ m_iterdata)
	{
	 string pepname;
	 set<string> prots;
	 double probvalue;
	 int chargestate;
	 
		// we have : multimap<string, pair< pair<double, int>, set<string> > > m_data;
		pepname     = m_iterdata->first;
		probvalue   = m_iterdata->second.first.first;
		chargestate = m_iterdata->second.first.second;
		prots       = m_iterdata->second.second;
		
		os << "e " << pepname << endl;
		os << "c " << chargestate << endl;
		
		set<string>::iterator setiter;
		for(setiter = prots.begin(); setiter != prots.end();++setiter)
		{
			os << "r " << *setiter << endl;
		}
		os << "p " << probvalue << endl;
		
	}
	
	os << "# end of automatically generated file" << endl;
	
}

void InputFileShrinker::read(istream & is)
{
	
	m_data.clear();
	m_chargeStates.clear();
	
  char instr;
  string pepName = "";
  string protName= "";
  double probvalue = -10;
  set<string> collectedProts;
  
  collectedProts.clear();
  
  int chargeState = 0;
  int state = 'e';


  while ( is >> instr )
   {
      //      cout << "Processing line: " << instr << endl;
      if ( instr == 'd' )
	  {
	   int charge;
	   double prior;
			is >> charge >> prior;
	  
		    //save charge state
			m_chargeStates[charge] = prior;
	  }
      else if ( instr == 'e' && (state == 'e' || state == 'p') )
	  {
	    if ( state == 'p' )
	    {
	      cerr << "Warning: no peptide score for peptide entry " << pepName << ", using last score (" << probvalue << ")" << endl;
	      
	      if ( pepName.length() <= 0 )
		  {
		    cerr << "Error: No previous peptide entry to use" << endl;
		    throw FormatException();
		  }
	    }

		// before we read the new peptide let's store the one we just read. 
		if(pepName.length() >0){
			m_data.insert( 
				pair<string, pair< pair<double, int>, set<string> > > (
					pepName, pair< pair<double, int>, set<string> > (
						pair<double, int>(probvalue, chargeState)
						, collectedProts
					)
			 ));
			 collectedProts.clear();
		}
	  
	    is >> pepName;

	  	state = 'c';
	  

	  }
      else if ( instr == 'c' && state == 'c' )
	  {
	    is >> chargeState;
	  
	    state = 'r';
	    
	  }
      else if ( instr == 'r' && ( state == 'c' || state == 'r' || state == 'p' ) )
	  {
		if(state == 'c'){  
			// charge state is skipped in this case
			chargeState = -1;
		}
	    is >> protName;
	    
	    collectedProts.insert(protName);
	    
	   state = 'p';
	  }
      else if ( instr == 'p' && state == 'p' )
	  {
	    is >> probvalue;
	    state = 'e';
	  }
      else if ( instr == '#' )
	  {
	    string garbage;
	    getline(is, garbage);
	  }
      else
	  {
	    cerr << "unexpected instruction " << instr << " in state " << state << endl;
	    string garbage;
	    getline(is, garbage);
	    cerr << "the input line was: " << instr << garbage << endl;
	    throw FormatException();
	  }
    }
    
    // finished
    cerr << "shrinker: read " << m_data.size() << " proteins" << endl;
}


