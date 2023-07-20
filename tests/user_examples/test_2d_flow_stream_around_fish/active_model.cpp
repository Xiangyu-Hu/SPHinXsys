#include "active_model.h"

#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
ActiveModelSolid::ActiveModelSolid(Real rho0, Real youngs_modulus, Real poisson_ratio)
    : SaintVenantKirchhoffSolid(rho0, youngs_modulus, poisson_ratio)
{
    material_type_name_ = "ActiveModelSolid";
}
//=============================================================================================//
void ActiveModelSolid::initializeLocalParameters(BaseParticles *base_particles)
{
    SaintVenantKirchhoffSolid::initializeLocalParameters(base_particles);
    base_particles->registerVariable(active_strain_, "ActiveStrain");
}
//=================================================================================================//
Matd ActiveModelSolid::StressPK1(Matd &F, size_t index_i)
{
    // Calculate the active deformation gradient
    Matd F0 = (2.0 * active_strain_[index_i] + Matd::Identity()).llt().matrixL();
    Matd F0_inv = F0.inverse();

    // Calculate the elastic deformation
    Matd F_e = F * F0_inv;
    Matd F0_star = F0.determinant() * F0_inv.transpose();

    // Calculate the elastic strain
    Matd E_e = 0.5 * (F.transpose() * F - Matd::Identity()) - active_strain_[index_i];
    return F_e * (lambda0_ * E_e.trace() * Matd::Identity() + 2.0 * G0_ * E_e) * F0_star;
}
//=================================================================================================//
} // namespace SPH
//=================================================================================================//
