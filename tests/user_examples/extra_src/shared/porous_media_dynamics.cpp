
#include "porous_media_dynamics.h"

#include <numeric>

namespace SPH
{
//=========================================================================================================//
namespace multi_species_continuum
{
//=================================================================================================//
BasePorousMediaRelaxation::
    BasePorousMediaRelaxation(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()),
    PorousMediaSolidDataInner(inner_relation), Vol_(particles_->Vol_),  
    pos_(particles_->pos_), vel_(particles_->vel_), 
    B_(particles_->B_), F_(particles_->F_), dF_dt_(particles_->dF_dt_) 
    {
        rho0_ = particles_->porous_solid_.ReferenceDensity();
        inv_rho0_ = 1.0 / rho0_;
        smoothing_length_ =   sph_body_.sph_adaptation_->ReferenceSmoothingLength();
    }
//=================================================================================================//
MomentumConstraint::MomentumConstraint(BodyPartByParticle &body_part)
	: BaseLocalDynamics<BodyPartByParticle>(body_part), PorousMediaSolidDataSimple(body_part.getSPHBody()),
	total_momentum_(particles_->total_momentum_) {}
//=================================================================================================//
PorousMediaStressRelaxationFirstHalf::
    PorousMediaStressRelaxationFirstHalf(BaseInnerRelation &body_inner_relation):
    BasePorousMediaRelaxation(body_inner_relation),
    Vol_update_(particles_->Vol_update_), fluid_saturation_(particles_->fluid_saturation_),
    total_mass_(particles_->total_mass_), fluid_mass_(particles_->fluid_mass_),
    dfluid_mass_dt_(particles_->dfluid_mass_dt_), total_momentum_(particles_->total_momentum_),
    dtotal_momentum_dt_(particles_->dtotal_momentum_dt_), 
    fluid_velocity_(particles_->fluid_velocity_),relative_fluid_flux_(particles_->relative_fluid_flux_),
    outer_fluid_velocity_relative_fluid_flux_(particles_->outer_fluid_velocity_relative_fluid_flux_),
    Stress_(particles_->Stress_)
    {
        numerical_dissipation_factor_ = 0.25;	

        fluid_initial_density = particles_->porous_solid_.getFulidInitialDensity();
        diffusivity_constant = particles_->porous_solid_.getDiffusivityConstant();
        water_pressure_constant= particles_->porous_solid_.getWaterPressureConstant();
    }
//=================================================================================================//
void PorousMediaStressRelaxationFirstHalf::initialization(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * 0.5;
    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
    Real J = F_[index_i].determinant();
    Matd inverse_F_T = F_[index_i].inverse().transpose();
    Matd  almansi_strain = 0.5 * (Matd::Identity() - (F_[index_i] * F_[index_i].transpose()).inverse());
    
    // first, update solid density and volume
    Vol_update_[index_i] = Vol_[index_i] * J;
        
    //obtain the Stress_PK one without J, J is added in total momentum
    //Note if anisotropic kernel is used, B correction cannot be applied here
    Stress_[index_i] = (particles_->porous_solid_.StressCauchy(almansi_strain, F_[index_i], index_i)
        - water_pressure_constant * (fluid_saturation_[index_i] - Eps) * Matd::Identity() ) * inverse_F_T;

    outer_fluid_velocity_relative_fluid_flux_[index_i] = fluid_velocity_[index_i] * relative_fluid_flux_[index_i].transpose() * inverse_F_T;
          
} 
//=================================================================================================//
void PorousMediaStressRelaxationFirstHalf::update(size_t index_i, Real dt)
{
    total_momentum_[index_i] += dtotal_momentum_dt_[index_i] * dt;
}
//=================================================================================================//
void PorousMediaStressRelaxationSecondHalf::initialization(size_t index_i, Real dt)
{
    // update solid velocity.
    vel_[index_i] = (total_momentum_[index_i] - relative_fluid_flux_[index_i]) * Vol_update_[index_i] / total_mass_[index_i];
    // update fluid velocity based on updated relative velocity, fluid density and solid velocity,
    fluid_velocity_[index_i] = vel_[index_i] - relative_fluid_flux_[index_i] / fluid_initial_density / (fluid_saturation_[index_i]+ TinyReal);
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
    fluid_saturation_[index_i] = fluid_mass_[index_i] / fluid_initial_density / Vol_update_[index_i];
}
//=================================================================================================// 
}
} // namespace SPH
