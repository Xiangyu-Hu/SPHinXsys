#include "base_material.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
void BaseMaterial::setLocalParameters(bool is_reload, BaseParticles *base_particles)
{
    if (!is_reload)
    {
        registerReloadLocalParameters(base_particles);
    }

    initializeLocalParameters(base_particles);
}
//=================================================================================================//
Fluid::Fluid(Real rho0, Real c0, Real mu)
    : BaseMaterial(rho0), c0_(c0), mu_(mu)
{
    material_type_name_ = "Fluid";
}
//=================================================================================================//
void Fluid::initializeLocalParameters(BaseParticles *base_particles)
{
    BaseMaterial::initializeLocalParameters(base_particles);
    base_particles->registerSharedVariable<Real>("Pressure", getPressure(rho0_));
}
//=================================================================================================//
StdLargeVec<Vecd> *Solid::AverageVelocity(BaseParticles *base_particles)
{
    return base_particles->getVariableByName<Vecd>("Velocity");
}
//=================================================================================================//
StdLargeVec<Vecd> *Solid::AverageForce(BaseParticles *base_particles)
{
    return base_particles->getVariableByName<Vecd>("Force");
}
//=================================================================================================//
} // namespace SPH
