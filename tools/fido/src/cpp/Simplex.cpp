#include "Simplex.h"

Simplex::Simplex(char * fname)
{
  ifstream fin(fname);

  load(fin);

  clear();
}

void Simplex::clear()
{
  buildTableau();

  basisColumns = (Array<int>) Set::FullSet(N, N+K-1);
  basisRows = (Array<int>) Set::FullSet(1, K);
}

void Simplex::buildTableau()
{
  N = f.size();
  K = A.numRows();
  
  mat = Matrix(K+1, N+K+1);

  Vector row = -Vector(f);
  row.resize(N+K+1);
  mat[0] = row;

  int k;
  for (k=0; k<K; k++)
    {
  // by negating, change the problem into A x <= b
      row = -A[k];

      row.resize(N+K+1);
      row.addElement(k+N, 1.0);

  // by negating, change the problem into A x <= b
      row.addElement(K+N, -b[k] );

      mat[k+1] = row;
      //      mat.add(row);
    }
}

void Simplex::pivot(int r, int c)
{
  int k;

  double cellValueRC = mat[r][c];

  mat[r] /= cellValueRC;

  for (k=0; k<=K; k++)
    {
      if ( k != r )
	{
	  double rowValue = mat[k][c];

	  if ( Vector::sparseChecker.isNonzero(rowValue) )
	    {
	      mat[k].addEqScaled( -rowValue, mat[r] );
	    }
	}
    }

  updateBasisIndices(r, c);
}

void Simplex::updateBasisIndices(int r, int c)
{
  // remove the column and row that was at the row
  Set st = basisRows == r;
  if ( st.size() != 1 )
    {
      cerr << "Error in basis indexing" << endl;
      exit(1);
    }
  
  basisRows.remove( st[0] );
  basisColumns.remove( st[0] );

  // add r and c to the rows and columns
  basisRows.add(r);
  basisColumns.add(c);
}

Vector Simplex::solve()
{
  int col;
  for ( iterations = 1; ; iterations++)
    {
      // removed for timing
      if ( iterations % 100 == 0 )
	{
	  //	  cout << "Objective " << mat[0][N+K] << endl;
	}

      Vector edges = mat[0][ Set::FullSet(0, N+K) ];

      if ( (edges < 0 ).isEmpty())
	{
	  // removed for timing
	  //	  cerr << "No negative element found in row 0" << endl;
	  break;
	}

      col = edges.argmin()[0];

      //      cout << "\tNeg element was " << mat[0].getElement(col) << endl;

      
      Vector column = mat.column(col);

      Set pos = column > 0.0;
      if ( pos.size() == 0 )
	{
	  // removed for timing
	  cerr << "No value in column found that is > 0" << endl;
	  break;
	}

      // find the pivot row offering greatest benefit

      //      cout << "Getting best column" << endl;
      Vector lastColumn = mat.column(N+K);
      int bestRow = pos[ (lastColumn / column)[pos].argmin() ][0];

      //      cout << "Pivoting " << bestRow << ", " << col << endl;
      pivot(bestRow, col);
    }

  Array<double> xStarArray(N, 0.0);
  for (int j=0; j<basisColumns.size(); j++)
    {
      if ( basisColumns[j] < N )
	{
	  // it is not a slack variable
	  xStarArray[ basisColumns[j] ] = mat[ basisRows[j] ][N+K];
	}
    }
  x = xStarArray;
  cout << names[ x > 0.0 ] << endl;

  //  cout << sol << endl;
  //  cout << "Objective value is : " << mat[0][N+K] << endl;

  // note: fixme. This should be the optimal x value
  return Vector();
}

