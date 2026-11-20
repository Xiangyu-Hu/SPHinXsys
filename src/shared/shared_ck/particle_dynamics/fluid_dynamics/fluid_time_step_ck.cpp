#include "fluid_time_step_ck.h"

#include "adaptation.h"
#include "viscosity.h"
#include "base_body.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
AdvectionTimeStepCK::
    AdvectionTimeStepCK(SPHBody &sph_body, Real U_ref, Real advectionCFL)
    : LocalDynamicsReduce<ReduceMax>(sph_body),
      h_min_(sph_body.getSPHAdaptation().MinimumSmoothingLength()),
      speed_ref_(U_ref), advectionCFL_(advectionCFL),
      dv_vel_(particles_->getVariableByName<Vecd>("Velocity")) {}
//=================================================================================================//
AdvectionTimeStepCK::FinishDynamics::FinishDynamics(AdvectionTimeStepCK &encloser)
    : encloser_(encloser) {}
//=================================================================================================//
Real AdvectionTimeStepCK::FinishDynamics::Result(Real reduced_value)
{
    // Live speed_ref_: AdvectionViscousTimeStepCK may refresh after μ(T) updates.
    return encloser_.advectionCFL_ * encloser_.h_min_ /
           (SMAX(std::sqrt(reduced_value), encloser_.speed_ref_) + TinyReal);
}
//=================================================================================================//
AdvectionViscousTimeStepCK::AdvectionViscousTimeStepCK(SPHBody &sph_body, Real U_ref, Real advectionCFL)
    : AdvectionTimeStepCK(sph_body, U_ref, advectionCFL), u_ref_(U_ref)
{
    refreshViscousSpeedReference();
}
//=================================================================================================//
void AdvectionViscousTimeStepCK::refreshViscousSpeedReference()
{
    auto &material = sph_body_->getMatterMaterial();
    Viscosity &viscosity = sph_body_->getMaterialProperty<Viscosity>();
    Real viscous_speed = viscosity.ReferenceViscosity() / material.ReferenceDensity() / h_min_;
    speed_ref_ = SMAX(viscous_speed, u_ref_);
}
//=================================================================================================//
AdvectionStepSetup::AdvectionStepSetup(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      dv_Vol_(particles_->getVariableByName<Real>("VolumetricMeasure")),
      dv_mass_(particles_->getVariableByName<Real>("Mass")),
      dv_rho_(particles_->getVariableByName<Real>("Density")),
      dv_dpos_(particles_->registerStateVariable<Vecd>("Displacement")) {}
//=================================================================================================//
UpdateParticlePosition::UpdateParticlePosition(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_dpos_(particles_->getVariableByName<Vecd>("Displacement")) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
