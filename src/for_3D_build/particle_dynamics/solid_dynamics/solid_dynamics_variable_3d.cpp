#include "base_body.h"
#include "solid_dynamics_variable.h"

#include <iterator>

namespace SPH
{
//=================================================================================================//
void VonMisesStrain::update(size_t index_i, Real dt)
{
    Mat3d F = F_[index_i];
    Mat3d green_lagrange_strain = 0.5 * (F.transpose() * F - Matd::Identity());

    Real epsilonxx = green_lagrange_strain(0, 0);
    Real epsilonyy = green_lagrange_strain(1, 1);
    Real epsilonzz = green_lagrange_strain(2, 2);
    Real epsilonxy = green_lagrange_strain(0, 1);
    Real epsilonxz = green_lagrange_strain(0, 2);
    Real epsilonyz = green_lagrange_strain(1, 2);

    derived_variable_[index_i] = sqrt((1.0 / 3.0) * (pow(epsilonxx - epsilonyy, 2) +
                                                     pow(epsilonyy - epsilonzz, 2) + pow(epsilonzz - epsilonxx, 2)) +
                                      2.0 * (pow(epsilonxy, 2) + pow(epsilonyz, 2) + pow(epsilonxz, 2)));
}
//=============================================================================================//
void VonMisesStrainDynamic::update(size_t index_i, Real dt)
{
    Mat3d F = F_[index_i];
    Mat3d green_lagrange_strain = 0.5 * (F.transpose() * F - Matd::Identity());

    Vec3d principal_strains = getPrincipalValuesFromMatrix(green_lagrange_strain);
    Real eps_1 = principal_strains[0];
    Real eps_2 = principal_strains[1];
    Real eps_3 = principal_strains[2];

    derived_variable_[index_i] = 1.0 / (1.0 + poisson_ratio_) *
                                 std::sqrt(0.5 * (pow(eps_1 - eps_2, 2) + pow(eps_2 - eps_3, 2) + pow(eps_3 - eps_1, 2)));
}
//=============================================================================================//
void VonMisesStress::update(size_t index_i, Real dt)
{
    Real J = rho0_ / rho_[index_i];
    Matd F = F_[index_i];
    Matd stress_PK1 = F * elastic_solid_.StressPK2(F, index_i);
    Matd sigma = (stress_PK1 * F.transpose()) / J;

    Real sigmaxx = sigma(0, 0);
    Real sigmayy = sigma(1, 1);
    Real sigmazz = sigma(2, 2);
    Real sigmaxy = sigma(0, 1);
    Real sigmaxz = sigma(0, 2);
    Real sigmayz = sigma(1, 2);

    derived_variable_[index_i] =
        sqrt(sigmaxx * sigmaxx + sigmayy * sigmayy + sigmazz * sigmazz -
             sigmaxx * sigmayy - sigmaxx * sigmazz - sigmayy * sigmazz +
             3.0 * (sigmaxy * sigmaxy + sigmaxz * sigmaxz + sigmayz * sigmayz));
}
//=============================================================================================//
void MidSurfaceVonMisesStress::update(size_t index_i, Real dt)
{
    Matd sigma = mid_surface_cauchy_stress_[index_i];

    Real sigmaxx = sigma(0, 0);
    Real sigmayy = sigma(1, 1);
    Real sigmazz = sigma(2, 2);
    Real sigmaxy = sigma(0, 1);
    Real sigmaxz = sigma(0, 2);
    Real sigmayz = sigma(1, 2);

    derived_variable_[index_i] =
        sqrt(sigmaxx * sigmaxx + sigmayy * sigmayy + sigmazz * sigmazz -
             sigmaxx * sigmayy - sigmaxx * sigmazz - sigmayy * sigmazz +
             3.0 * (sigmaxy * sigmaxy + sigmaxz * sigmaxz + sigmayz * sigmayz));
}
//=================================================================================================//
} // namespace SPH
