#include "base_material.h"

#include "base_particles.hpp"
#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
BaseMaterial::~BaseMaterial() = default;
//=================================================================================================//
void BaseMaterial::setLocalParameters(SPHSystem &sph_system, BaseParticles *base_particles)
{
    if (sph_system.ReloadParticles())
    {
        registerLocalParametersFromReload(base_particles);
    }
    else
    {
        registerLocalParametersToReload(base_particles);
    }

    initializeLocalParameters(base_particles);
}
//=================================================================================================//
std::string BaseMaterial::Name()
{
    if (material_name_.empty())
    {
        return material_type_name_;
    }
    return material_name_;
}
//=================================================================================================//
MatterMaterial::MatterMaterial() : BaseMaterial()
{
    material_type_name_ = "MatterMaterial";
}
//=================================================================================================//
void MatterMaterial::initializeLocalParameters(BaseParticles *base_particles)
{
    dv_rho_ = base_particles->registerStateVariable<Real>("Density", ReferenceDensity());
    dv_mass_ = base_particles->registerStateVariable<Real>(
        "Mass", [&](UnsignedInt i) -> Real
        { return ReferenceDensity() * base_particles->ParticleVolume(i); });
}
//=================================================================================================//
Fluid::Fluid() : MatterMaterial()
{
    material_type_name_ = "Fluid";
}
//=================================================================================================//
SolidContact::SolidContact(Real rho0, Real contact_stiffness, Real contact_friction)
    : rho0_copy_(rho0), contact_stiffness_(contact_stiffness), contact_friction_(contact_friction) {}
//=================================================================================================//
Solid::Solid(Real rho0, Real contact_stiffness, Real contact_friction)
    : MatterMaterial(), SolidContact(rho0, contact_stiffness, contact_friction), rho0_(rho0)
{
    material_type_name_ = "Solid";
}
//=================================================================================================//
Vecd *Solid::AverageVelocity(BaseParticles *base_particles)
{
    return base_particles->registerStateVariableData<Vecd>("Velocity");
}
//=================================================================================================//
Vecd *Solid::AverageAcceleration(BaseParticles *base_particles)
{
    return base_particles->registerStateVariableData<Vecd>("Acceleration");
}
//=================================================================================================//
DiscreteVariable<Vecd> *Solid::AverageVelocityVariable(BaseParticles *base_particles)
{
    return base_particles->registerStateVariable<Vecd>("Velocity");
}
//=================================================================================================//
DiscreteVariable<Vecd> *Solid::AverageAccelerationVariable(BaseParticles *base_particles)
{
    return base_particles->registerStateVariable<Vecd>("Acceleration");
}
//=================================================================================================//
} // namespace SPH
