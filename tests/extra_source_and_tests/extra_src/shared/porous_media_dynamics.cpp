
#include "porous_media_dynamics.h"

#include <numeric>

namespace SPH
{
namespace multi_species_continuum
{
//=================================================================================================//
GetSaturationTimeStepSize::GetSaturationTimeStepSize(SPHBody &sph_body)
    : LocalDynamicsReduce<ReduceMin>(sph_body),
      porous_solid_(DynamicCast<PorousMediaSolid>(this, particles_->getBaseMaterial())),
      smoothing_length_(sph_body.sph_adaptation_->ReferenceSmoothingLength()) {}
//=================================================================================================//
BasePorousMediaRelaxation::BasePorousMediaRelaxation(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()),
      DataDelegateInner(inner_relation),
      porous_solid_(DynamicCast<PorousMediaSolid>(this, particles_->getBaseMaterial())),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      vel_(particles_->registerSharedVariable<Vecd>("Velocity")),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
      F_(particles_->registerSharedVariable<Matd>("DeformationGradient", IdentityMatrix<Matd>::value)),
      dF_dt_(particles_->registerSharedVariable<Matd>("DeformationRate"))
{
    rho0_ = porous_solid_.ReferenceDensity();
    inv_rho0_ = 1.0 / rho0_;
    smoothing_length_ = sph_body_.sph_adaptation_->ReferenceSmoothingLength();
}
//=================================================================================================//
MomentumConstraint::MomentumConstraint(BodyPartByParticle &body_part)
    : BaseLocalDynamics<BodyPartByParticle>(body_part),
      total_momentum_(particles_->getVariableDataByName<Vecd>("TotalMomentum")) {}
//=================================================================================================//
PorousMediaStressRelaxationFirstHalf::
    PorousMediaStressRelaxationFirstHalf(BaseInnerRelation &body_inner_relation)
    : BasePorousMediaRelaxation(body_inner_relation),
      Vol_update_(particles_->registerSharedVariable<Real>("UpdateVolume")),
      fluid_saturation_(particles_->registerSharedVariable<Real>("FluidSaturation")),
      total_mass_(particles_->registerSharedVariable<Real>("TotalMass")),
      fluid_mass_(particles_->registerSharedVariable<Real>("FluidMass")),
      dfluid_mass_dt_(particles_->registerSharedVariable<Real>("FluidMassIncrement")),
      total_momentum_(particles_->registerSharedVariable<Vecd>("TotalMomentum")),
      force_(particles_->registerSharedVariable<Vecd>("Force")),
      force_prior_(particles_->registerSharedVariable<Vecd>("ForcePrior")),
      fluid_velocity_(particles_->registerSharedVariable<Vecd>("FluidVelocity")),
      relative_fluid_flux_(particles_->registerSharedVariable<Vecd>("RelativeFluidFlux")),
      outer_fluid_velocity_relative_fluid_flux_(particles_->registerSharedVariable<Matd>("OuterFluidVelocityRelativeFluidFlux")),
      Stress_(particles_->registerSharedVariable<Matd>("Stress")),
      diffusivity_constant_(porous_solid_.getDiffusivityConstant()),
      fluid_initial_density_(porous_solid_.getFluidInitialDensity()),
      water_pressure_constant_(porous_solid_.getWaterPressureConstant()) {}
//=================================================================================================//
void PorousMediaStressRelaxationFirstHalf::initialization(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * 0.5;
    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
    Real J = F_[index_i].determinant();
    Matd inverse_F_T = F_[index_i].inverse().transpose();
    Matd almansi_strain = 0.5 * (Matd::Identity() - (F_[index_i] * F_[index_i].transpose()).inverse());

    // first, update solid density and volume
    Vol_update_[index_i] = Vol_[index_i] * J;

    // obtain the Stress_PK one without J, J is added in total momentum
    // Note if anisotropic kernel is used, B correction cannot be applied here
    Stress_[index_i] = (porous_solid_.StressCauchy(almansi_strain, index_i) - water_pressure_constant_ * (fluid_saturation_[index_i] - Eps) * Matd::Identity()) * inverse_F_T;

    outer_fluid_velocity_relative_fluid_flux_[index_i] = fluid_velocity_[index_i] * relative_fluid_flux_[index_i].transpose() * inverse_F_T;
}
//=================================================================================================//
void PorousMediaStressRelaxationFirstHalf::update(size_t index_i, Real dt)
{
    total_momentum_[index_i] += (force_prior_[index_i] + force_[index_i]) * dt;
}
//=================================================================================================//
void PorousMediaStressRelaxationSecondHalf::initialization(size_t index_i, Real dt)
{
    // update solid velocity.
    vel_[index_i] = (total_momentum_[index_i] - relative_fluid_flux_[index_i]) * Vol_update_[index_i] / total_mass_[index_i];
    // update fluid velocity based on updated relative velocity, fluid density and solid velocity,
    fluid_velocity_[index_i] = vel_[index_i] - relative_fluid_flux_[index_i] / fluid_initial_density_ / (fluid_saturation_[index_i] + TinyReal);
    pos_[index_i] += vel_[index_i] * dt * 0.5;
}
//=================================================================================================//
void PorousMediaStressRelaxationSecondHalf::update(size_t index_i, Real dt)
{
    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
}
//=================================================================================================//
void SaturationRelaxationInPorousMedia::initialization(size_t index_i, Real Dt) {}
//=================================================================================================//
void SaturationRelaxationInPorousMedia::update(size_t index_i, Real Dt)
{
    fluid_mass_[index_i] += dfluid_mass_dt_[index_i] * Dt;
    // update total mass
    total_mass_[index_i] = rho0_ * Vol_[index_i] + fluid_mass_[index_i];
    //  update fluid saturation
    fluid_saturation_[index_i] = fluid_mass_[index_i] / fluid_initial_density_ / Vol_update_[index_i];
}
//=================================================================================================//
} // namespace multi_species_continuum
} // namespace SPH
