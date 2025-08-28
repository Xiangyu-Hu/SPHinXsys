#include "base_particles.hpp"
#include "elastic_solid.h"

namespace SPH
{
//=================================================================================================//
void OrthotropicSolid::CalculateAllMu()
{

    Mu_[0] = 1 / G_[0] + 1 / G_[2] - 1 / G_[1];
    Mu_[1] = 1 / G_[1] + 1 / G_[0] - 1 / G_[2];
    Mu_[2] = 1 / G_[2] + 1 / G_[1] - 1 / G_[0];
}
//=================================================================================================//
void OrthotropicSolid::CalculateAllLambda()
{
    // first we calculate the upper left part, a 3x3 matrix of the full compliance matrix
    Matd Compliance = Matd::Zero();
    Compliance.col(0) = Vecd(1 / E_[0], -poisson_[0] / E_[1], -poisson_[1] / E_[2]);
    Compliance.col(1) = Vecd(-poisson_[0] / E_[0], 1 / E_[1], -poisson_[2] / E_[1]);
    Compliance.col(2) = Vecd(-poisson_[1] / E_[0], -poisson_[2] / E_[1], 1 / E_[2]);
    // we calculate the inverse of the Compliance matrix, and calculate the lambdas elementwise
    Matd Compliance_inv = Compliance.inverse();
    // Lambda_ is a 3x3 matrix
    Lambda_(0, 0) = Compliance_inv(0, 0) - 2 * Mu_[0];
    Lambda_(1, 1) = Compliance_inv(1, 1) - 2 * Mu_[1];
    Lambda_(2, 2) = Compliance_inv(2, 2) - 2 * Mu_[2];
    Lambda_(0, 1) = Compliance_inv(0, 1);
    Lambda_(0, 2) = Compliance_inv(0, 2);
    Lambda_(1, 2) = Compliance_inv(1, 2);
    // the matrix is symmetric
    Lambda_(1, 0) = Lambda_(0, 1);
    Lambda_(2, 0) = Lambda_(0, 2);
    Lambda_(2, 1) = Lambda_(1, 2);
}
//=================================================================================================//
} // namespace SPH