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
void BaseMaterial::initializeLocalParameters(BaseParticles *base_particles)
{
    base_particles->registerSharedVariable<Vecd>("Velocity");
    base_particles->registerSharedVariable<Vecd>("Force");
    base_particles->registerSharedVariable<Vecd>("PriorForce");
    base_particles->registerSharedVariable("Density", rho0_);
    base_particles->registerSharedVariable<Real>("Mass",
                                                 [&](size_t i) -> Real
                                                 { return rho0_ * base_particles->ParticleVolume(i); });
    //----------------------------------------------------------------------
    //		add particle data for output
    //----------------------------------------------------------------------
    base_particles->addVariableToWrite<Vecd>("Velocity");
    //----------------------------------------------------------------------
    //		add particle data for restart
    //----------------------------------------------------------------------
    base_particles->addVariableToRestart<Vecd>("Velocity");
    base_particles->addVariableToRestart<Vecd>("Force");
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
    base_particles->registerSharedVariable("Pressure", getPressure(rho0_));
}
//=================================================================================================//
} // namespace SPH
