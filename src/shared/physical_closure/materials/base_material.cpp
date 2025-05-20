#include "base_material.h"

#include "base_particles.hpp"
#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
void BaseMaterial::setLocalParameters(SPHSystem &sph_system, BaseParticles *base_particles)
{
    if (sph_system.ReloadParticles())
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
Fluid::Fluid(Real rho0, Real c0) : BaseMaterial(rho0), c0_(c0)
{
    material_type_name_ = "Fluid";
}
//=================================================================================================//
Fluid::EosKernel::EosKernel(Fluid &encloser) : c0_(encloser.c0_), rho0_(encloser.rho0_) {}
//=================================================================================================//
Vecd *Solid::AverageVelocity(BaseParticles *base_particles)
{
    return base_particles->registerStateVariable<Vecd>("Velocity");
}
//=================================================================================================//
Vecd *Solid::AverageAcceleration(BaseParticles *base_particles)
{
    return base_particles->registerStateVariable<Vecd>("Acceleration");
}
//=================================================================================================//
DiscreteVariable<Vecd> *Solid::AverageVelocityVariable(BaseParticles *base_particles)
{
    return base_particles->registerStateVariableOnly<Vecd>("Velocity");
}
//=================================================================================================//
DiscreteVariable<Vecd> *Solid::AverageAccelerationVariable(BaseParticles *base_particles)
{
    return base_particles->registerStateVariableOnly<Vecd>("Acceleration");
}
//=================================================================================================//
} // namespace SPH
