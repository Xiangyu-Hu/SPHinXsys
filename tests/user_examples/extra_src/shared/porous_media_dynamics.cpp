
#include "porous_media_dynamics.h"

#include <numeric>

namespace SPH
{
//=========================================================================================================//
namespace solid_dynamics
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
void PorousMediaStressRelaxationFirstHalf::interaction(size_t index_i, Real dt)
{
    Vecd total_momentum_increment = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];

    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real r_ij = inner_neighborhood.r_ij_[n];
        Vecd gradw_ijV_j_ = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
 
        Real dim_r_ij_1 = Dimensions / r_ij;
        Vecd pos_jump = pos_[index_i] - pos_[index_j];
        Vecd vel_jump = vel_[index_i] - vel_[index_j];
        Real strain_rate =  pos_jump.dot(vel_jump) * dim_r_ij_1 * dim_r_ij_1;
        Real weight = inner_neighborhood.W_ij_[n] * inv_W0_;

        Matd numerical_stress_ij = 0.5 * (F_[index_i] + F_[index_j]) 
            * particles_->porous_solid_.PairNumericalDamping(strain_rate,  smoothing_length_);
    
        //three parts for the momentum increment
        total_momentum_increment += (Stress_[index_i] + Stress_[index_j])* gradw_ijV_j_
                                    +  numerical_dissipation_factor_  * numerical_stress_ij * weight * gradw_ijV_j_ ;
    
        total_momentum_increment -= (outer_fluid_velocity_relative_fluid_flux_[index_i] + outer_fluid_velocity_relative_fluid_flux_[index_j])
                                     *gradw_ijV_j_;
    }

    dtotal_momentum_dt_[index_i] = total_momentum_increment;
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
void PorousMediaStressRelaxationSecondHalf::interaction(size_t index_i, Real dt)
{
    Matd deformation_gradient_change_rate = Matd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd gradw_ijV_j_ = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];

        deformation_gradient_change_rate -=
               (vel_[index_i] - vel_[index_j]) * gradw_ijV_j_.transpose();
    }
    dF_dt_[index_i] = deformation_gradient_change_rate * B_[index_i];
}
//=================================================================================================//
void PorousMediaStressRelaxationSecondHalf::update(size_t index_i, Real dt)
{
    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
}
//=================================================================================================//
void SaturationRelaxationInPorousMedia::initialization(size_t index_i, Real Dt) {}
//=================================================================================================//
void SaturationRelaxationInPorousMedia::interaction(size_t index_i, Real Dt) 
{
    Vecd fluid_saturation_gradient  = Vecd::Zero();
    Real relative_fluid_flux_divergence  = 0.0;	
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];

    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real r_ij = inner_neighborhood.r_ij_[n];
        Real dw_ijV_j_ = inner_neighborhood.dW_ijV_j_[n] ;

        Vecd e_ij = inner_neighborhood.e_ij_[n];
        fluid_saturation_gradient -=  (fluid_saturation_[index_i] - fluid_saturation_[index_j])
                                        * e_ij* dw_ijV_j_;
                
        
        relative_fluid_flux_divergence +=  1.0 / 2.0 * (fluid_saturation_[index_i] * fluid_saturation_[index_i] - fluid_saturation_[index_j] * fluid_saturation_[index_j]) 
                                                / (r_ij + TinyReal) *  dw_ijV_j_;
                                        
    }
    //then we update relative velocity based on the updated fluid density
    relative_fluid_flux_[index_i] = -diffusivity_constant * fluid_initial_density * fluid_saturation_[index_i] * fluid_saturation_gradient;
    
    dfluid_mass_dt_[index_i] = diffusivity_constant * Vol_update_[index_i] * fluid_initial_density
                    * relative_fluid_flux_divergence;
}		
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
} // namespace solid_dynamics
} // namespace SPH
