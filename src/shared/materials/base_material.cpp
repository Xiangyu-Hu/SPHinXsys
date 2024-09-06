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
    return base_particles->registerDiscreteVariablesOnly<Vecd>("Velocity");
}
//=================================================================================================//
DiscreteVariable<Vecd> *Solid::AverageAccelerationVariable(BaseParticles *base_particles)
{
    return base_particles->registerDiscreteVariablesOnly<Vecd>("Acceleration");
}
//=================================================================================================//
} // namespace SPH
