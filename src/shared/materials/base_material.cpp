#include "base_material.h"

#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
void BaseMaterial::setLocalParameters(bool is_reload, BaseParticles *base_particles)
{
    if (is_reload)
    {
        registerLocalParametersFromReload(base_particles);
    }
    else
    {
        registerLocalParameters(base_particles);
    }

    initializeLocalParameters(base_particles);
}
//=================================================================================================//
Fluid::Fluid(Real rho0, Real c0, Real mu) : BaseMaterial(rho0), c0_(c0), mu_(mu)
{
    material_type_name_ = "Fluid";
}
//=================================================================================================//
StdLargeVec<Vecd> *Solid::AverageVelocity(BaseParticles *base_particles)
{
    return base_particles->registerSharedVariable<Vecd>("Velocity");
}
//=================================================================================================//
StdLargeVec<Vecd> *Solid::AverageAcceleration(BaseParticles *base_particles)
{
    return base_particles->registerSharedVariable<Vecd>("Acceleration");
}
//=================================================================================================//
} // namespace SPH
