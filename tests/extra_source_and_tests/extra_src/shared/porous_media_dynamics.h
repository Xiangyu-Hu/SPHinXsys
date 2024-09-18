/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	POROUS_ELASTIC_DYNAMICS_H.h
 * @brief 	Here, we define the algorithm classes for elastic solid dynamics.
 * @details 	We consider here a weakly compressible solids.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef POROUS_ELASTIC_DYNAMICS_H
#define POROUS_ELASTIC_DYNAMICS_H

#include "all_particle_dynamics.h"
#include "elastic_dynamics.h"
#include "porous_media_solid.h"

namespace SPH
{
namespace multi_species_continuum
{
/**
 * @class GetSaturationTimeStepSize
 * @brief Computing the time step size based on diffusion coefficient and particle smoothing length
 */
class GetSaturationTimeStepSize
    : public LocalDynamicsReduce<ReduceMin>
{
  protected:
    PorousMediaSolid &porous_solid_;
    Real saturation_time_step_;
    Real smoothing_length_;

  public:
    explicit GetSaturationTimeStepSize(SPHBody &sph_body);
    virtual ~GetSaturationTimeStepSize(){};

    Real reduce(size_t index_i, Real dt = 0.0)
    {
        return 0.5 * smoothing_length_ * smoothing_length_ /
               porous_solid_.getDiffusivityConstant() / (Real)Dimensions;
    };
};

/**@class MomentumConstraint
 * @brief MomentumConstraint with zero momentum.
 */
class MomentumConstraint : public BaseLocalDynamics<BodyPartByParticle>
{
  public:
    explicit MomentumConstraint(BodyPartByParticle &body_part);
    virtual ~MomentumConstraint(){};

    void update(size_t index_i, Real dt = 0.0) { total_momentum_[index_i] = Vecd::Zero(); };

  protected:
    Vecd *total_momentum_;
};

/**
 * @class BasePorousStressRelaxation
 * @brief base class for porous media relaxation
 */
class BasePorousMediaRelaxation : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit BasePorousMediaRelaxation(BaseInnerRelation &inner_relation);
    virtual ~BasePorousMediaRelaxation(){};

  protected:
    PorousMediaSolid &porous_solid_;
    Real *Vol_;
    Vecd *pos_, *vel_;
    Matd *B_, *F_, *dF_dt_;
    Real rho0_, inv_rho0_;
    Real smoothing_length_;
};

/**
 * @class PorousMediaStressRelaxationFirstHalf
 * @brief computing Porous Media stress relaxation process by verlet time stepping
 * This is the first step
 */
class PorousMediaStressRelaxationFirstHalf
    : public BasePorousMediaRelaxation
{
  public:
    PorousMediaStressRelaxationFirstHalf(BaseInnerRelation &body_inner_relation);
    virtual ~PorousMediaStressRelaxationFirstHalf(){};

  protected:
    Real *Vol_update_, *fluid_saturation_, *total_mass_, *fluid_mass_, *dfluid_mass_dt_;
    Vecd *total_momentum_, *force_, *force_prior_, *fluid_velocity_, *relative_fluid_flux_;
    Matd *outer_fluid_velocity_relative_fluid_flux_, *Stress_;

    Real diffusivity_constant_, fluid_initial_density_, water_pressure_constant_;

    const Real one_over_dimensions_ = 1.0 / (Real)Dimensions;
    Real numerical_dissipation_factor_ = 0.25;
    Real inv_W0_ = 1.0 / sph_body_.sph_adaptation_->getKernel()->W0(ZeroVecd);

    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd total_momentum_increment = Vecd::Zero();
        Neighborhood &inner_neighborhood = inner_configuration_[index_i];

        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Real r_ij = inner_neighborhood.r_ij_[n];
            Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];

            Real dim_r_ij_1 = Dimensions / r_ij;
            Vecd pos_jump = pos_[index_i] - pos_[index_j];
            Vecd vel_jump = vel_[index_i] - vel_[index_j];
            Real strain_rate = pos_jump.dot(vel_jump) * dim_r_ij_1 * dim_r_ij_1;
            Real weight = inner_neighborhood.W_ij_[n] * inv_W0_;

            Matd numerical_stress_ij = 0.5 * (F_[index_i] + F_[index_j]) *
                                       porous_solid_.PairNumericalDamping(strain_rate, smoothing_length_);

            // three parts for the momentum increment
            total_momentum_increment += (Stress_[index_i] + Stress_[index_j] + numerical_dissipation_factor_ * numerical_stress_ij * weight -
                                         outer_fluid_velocity_relative_fluid_flux_[index_i] - outer_fluid_velocity_relative_fluid_flux_[index_j]) *
                                        gradW_ijV_j;
        }

        force_[index_i] = total_momentum_increment;
    };
    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class PorousMediaStressRelaxationSecondHalf
 * @brief computing Porous Media stress relaxation process by verlet time stepping
 * This is the second step
 */
class PorousMediaStressRelaxationSecondHalf
    : public PorousMediaStressRelaxationFirstHalf
{
  public:
    PorousMediaStressRelaxationSecondHalf(BaseInnerRelation &body_inner_relation)
        : PorousMediaStressRelaxationFirstHalf(body_inner_relation){};
    virtual ~PorousMediaStressRelaxationSecondHalf(){};

  protected:
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0)
    {
        Matd deformation_gradient_change_rate = Matd::Zero();
        Neighborhood &inner_neighborhood = inner_configuration_[index_i];

        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];

            deformation_gradient_change_rate -=
                (vel_[index_i] - vel_[index_j]) * gradW_ijV_j.transpose();
        }
        dF_dt_[index_i] = deformation_gradient_change_rate * B_[index_i];
    };
    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class PorousMediaSaturationDynamicsInitialCondition
 * @brief Set initial condition  in porous media.
 */
class PorousMediaSaturationDynamicsInitialCondition : public BaseLocalDynamics<BodyPartByParticle>
{
  public:
    PorousMediaSaturationDynamicsInitialCondition(BodyPartByParticle &body_part)
        : BaseLocalDynamics<BodyPartByParticle>(body_part),
          fluid_mass_(particles_->getVariableDataByName<Real>("FluidMass")),
          fluid_saturation_(particles_->getVariableDataByName<Real>("FluidSaturation")),
          total_mass_(particles_->getVariableDataByName<Real>("TotalMass")),
          rho_n_(particles_->getVariableDataByName<Real>("Density")),
          Vol_update_(particles_->getVariableDataByName<Real>("UpdateVolume")),
          pos_(particles_->getVariableDataByName<Vecd>("Position")){};

    virtual ~PorousMediaSaturationDynamicsInitialCondition(){};

  protected:
    Real *fluid_mass_, *fluid_saturation_, *total_mass_, *rho_n_, *Vol_update_;
    Vecd *pos_;
};

/**
 * @class SaturationRelaxationInPorousMedia
 * @brief computing saturation relaxation process in porous media
 */
class SaturationRelaxationInPorousMedia
    : public PorousMediaStressRelaxationFirstHalf
{
  public:
    SaturationRelaxationInPorousMedia(BaseInnerRelation &body_inner_relation)
        : PorousMediaStressRelaxationFirstHalf(body_inner_relation){};
    virtual ~SaturationRelaxationInPorousMedia(){};

  protected:
    void initialization(size_t index_i, Real Dt = 0.0);
    void interaction(size_t index_i, Real Dt = 0.0)
    {
        Vecd fluid_saturation_gradient = Vecd::Zero();
        Real relative_fluid_flux_divergence = 0.0;
        Neighborhood &inner_neighborhood = inner_configuration_[index_i];

        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Real r_ij = inner_neighborhood.r_ij_[n];
            Real dw_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];

            Vecd e_ij = inner_neighborhood.e_ij_[n];
            fluid_saturation_gradient -= (fluid_saturation_[index_i] - fluid_saturation_[index_j]) * e_ij * dw_ijV_j;

            relative_fluid_flux_divergence += 1.0 / 2.0 * (fluid_saturation_[index_i] * fluid_saturation_[index_i] - fluid_saturation_[index_j] * fluid_saturation_[index_j]) / (r_ij + TinyReal) * dw_ijV_j;
        }
        // then we update relative velocity based on the updated fluid density
        relative_fluid_flux_[index_i] = -diffusivity_constant_ * fluid_initial_density_ * fluid_saturation_[index_i] * fluid_saturation_gradient;

        dfluid_mass_dt_[index_i] = diffusivity_constant_ * Vol_update_[index_i] * fluid_initial_density_ * relative_fluid_flux_divergence;
    };
    void update(size_t index_i, Real Dt = 0.0);
};

} // namespace multi_species_continuum
} // namespace SPH
#endif // POROUS_ELASTIC_DYNAMICS_H
