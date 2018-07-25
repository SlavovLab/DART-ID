#include "EmissionModel.h"

EmissionModel::EmissionModel(double a, double b, double g)
{
  alpha = a;
  beta = b;
  gamma = g;
}

double EmissionModel::connectedEmission() const
{
}

