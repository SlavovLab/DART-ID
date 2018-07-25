#include "GridEmissionModel.h"

double GridEmissionModel::shared_alpha_min, GridEmissionModel::shared_alpha_max, GridEmissionModel::shared_beta_min, GridEmissionModel::shared_beta_max, GridEmissionModel::shared_resolution, GridEmissionModel::shared_gamma;

void GridEmissionModel::initialize(double res, double g, double a_min, double a_max, double b_min, double b_max)
{
  GridEmissionModel::shared_resolution = res;
  GridEmissionModel::shared_gamma = g;

  GridEmissionModel::shared_alpha_min = a_min;
  GridEmissionModel::shared_alpha_max = a_max;
  GridEmissionModel::shared_beta_min = b_min;
  GridEmissionModel::shared_beta_max = b_max;


  //  cout << "resolution: " << resolution << " ; \t" << alpha_min << ", " << alpha_max << "\t" << beta_min << ", " << beta_max << endl;
}

pair<int, int> GridEmissionModel::matrixIndices() const
{
  return pair<int,int>(alphaIndex, betaIndex);
}
  
pair<int, int> GridEmissionModel::matrixRange() const
{
  int alphaMaxIndex, betaMaxIndex;
    
  alphaMaxIndex = int( (alpha_max - alpha_min + resolution) / resolution );
  betaMaxIndex = int( (beta_max - beta_min + resolution) / resolution );

  return pair<int,int>(alphaMaxIndex, betaMaxIndex);
}

bool GridEmissionModel::advance()
{
  /***
  count++;

  alphaIndex++;

  alpha += resolution;

  if ( alpha >= ( alpha_max ) - 1e-6 )
    {
      return false;
    }

  return true;
  ***/

  count++;

  alphaIndex++;

  alpha += resolution;

  //  cout << "alpha=" << alpha << " alpha_max=" << alpha_max << endl;
  if ( alpha >= ( alpha_max ) - 1e-6 )
    {
      //      cout << "alpha has passed/equaled alpha_max" << endl;

      betaIndex++;
      beta += resolution;

      // triangular prior: alpha >= beta
      //      alpha = beta;
      //      alphaIndex = betaIndex;

      // rectangular prior
      alphaIndex = 0;
      alpha = alpha_min;
    }

  if ( beta != beta_min && beta  >= ( beta_max ) - 1e-6 )
  //  if ( beta  >= ( beta_max ) - 1e-6 )
    {
      return false;
    }

  return true;
}
