#include "base_material.h"
#include "base_particles.hpp"

namespace SPH
{
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
} // namespace SPH
