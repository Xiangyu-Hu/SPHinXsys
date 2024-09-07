/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	k-epsilon_turbulent_model.hpp
 * @brief 	k-epsilon_turbulent_model.hpp
 * @author	 Xiangyu Hu
 */
#pragma once

#include "k-epsilon_turbulent_model.h"

namespace SPH
{
//=====================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
BaseTurtbulentModel<Base, DataDelegationType>::BaseTurtbulentModel(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      turbu_k_(*(this->particles_->template registerSharedVariable<Real>("TurbulenceKineticEnergy"))),
      turbu_epsilon_(*(this->particles_->template registerSharedVariable<Real>("TurbulentDissipation"))),
      turbu_mu_(*(this->particles_->template registerSharedVariable<Real>("TurbulentViscosity"))),
      turbu_strain_rate_(*(this->particles_->template registerSharedVariable<Matd>("TurbulentStrainRate"))),
      mu_(DynamicCast<Fluid>(this, this->particles_->getBaseMaterial()).ReferenceViscosity()),
      smoothing_length_(this->sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
      particle_spacing_min_(base_relation.getSPHBody().sph_adaptation_->MinimumSpacing()),
      rho_(*this->particles_->template getVariableDataByName<Real>("Density")),
      Vol_(*this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      vel_(*this->particles_->template getVariableDataByName<Vecd>("Velocity")),
      dimension_(2) {}
//A temporarily treatment for dimention
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
TKEnergyForce<Base, DataDelegationType>::
    TKEnergyForce(BaseRelationType &base_relation)
    : BaseTurtbulentModel<Base, DataDelegationType>(base_relation),
      force_(*this->particles_->template registerSharedVariable<Vecd>("Force")),
      mass_(*this->particles_->template getVariableDataByName<Real>("Mass")),
      indicator_(*this->particles_->template getVariableDataByName<int>("Indicator")),
      pos_(*this->particles_->template getVariableDataByName<Vecd>("Position")),
      turbu_k_(*this->particles_->template getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      test_k_grad_rslt_(*this->particles_->template registerSharedVariable<Vecd>("TkeGradResult")) {}
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
GetVelocityGradient<DataDelegationType>::
    GetVelocityGradient(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      Vol_(*this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      vel_(*this->particles_->template getVariableDataByName<Vecd>("Velocity")),
      pos_(*this->particles_->template getVariableDataByName<Vecd>("Position")),
      is_near_wall_P1_(*this->particles_->template getVariableDataByName<int>("IsNearWallP1")),
      is_near_wall_P2_(*this->particles_->template getVariableDataByName<int>("IsNearWallP2")),
      velocity_gradient_(*(this->particles_->template registerSharedVariable<Matd>("TurbulentVelocityGradient"))),
      velocity_gradient_wall(*(this->particles_->template registerSharedVariable<Matd>("Velocity_Gradient_Wall"))) {}
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
TurbuViscousForce<DataDelegationType>::TurbuViscousForce(BaseRelationType &base_relation)
    : ViscousForce<DataDelegationType>(base_relation),
      turbu_k_(*this->particles_->template getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_mu_(*this->particles_->template getVariableDataByName<Real>("TurbulentViscosity")),
      wall_Y_plus_(*this->particles_->template getVariableDataByName<Real>("WallYplus")),
      wall_Y_star_(*this->particles_->template getVariableDataByName<Real>("WallYstar")),
      velo_friction_(*this->particles_->template getVariableDataByName<Vecd>("FrictionVelocity")),
      y_p_(*this->particles_->template getVariableDataByName<Real>("Y_P")),
      is_near_wall_P2_(*this->particles_->template getVariableDataByName<int>("IsNearWallP2")),
      molecular_viscosity_(DynamicCast<Fluid>(this, this->particles_->getBaseMaterial()).ReferenceViscosity()),
      c0_(DynamicCast<Fluid>(this, this->particles_->getBaseMaterial()).ReferenceSoundSpeed()) {}
//=================================================================================================//
template <class DataDelegationType, class ParticleScope>
template <class BaseRelationType>
ExtraTransportForce<Base, DataDelegationType, ParticleScope>::
    ExtraTransportForce(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      rho_(*this->particles_->template getVariableDataByName<Real>("Density")), vel_(*this->particles_->template getVariableDataByName<Vecd>("Velocity")),
      zero_gradient_residue_(*this->particles_->template getVariableDataByName<Vecd>("ZeroGradientResidue")),
      extra_transport_stress_(*(this->particles_->template registerSharedVariable<Matd>("ExtraTransportStress"))),
      extra_transport_vel_(*(this->particles_->template registerSharedVariable<Vecd>("ExtraTransportVelocity"))),
      within_scope_(this->particles_)
{
    static_assert(std::is_base_of<WithinScope, ParticleScope>::value,
                  "WithinScope is not the base of ParticleScope!");
}

//=================================================================================================//
template <class LimiterType, typename... CommonControlTypes>
ExtraTransportForce<Inner<LimiterType>, CommonControlTypes...>::ExtraTransportForce(BaseInnerRelation &inner_relation)
    : ExtraTransportForce<Base, DataDelegateInner, CommonControlTypes...>(inner_relation),
      Vol_(*this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      h_ref_(this->sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
      extra_transport_stress_(*this->particles_->template getVariableDataByName<Matd>("ExtraTransportStress")),
      extra_transport_vel_(*this->particles_->template getVariableDataByName<Vecd>("ExtraTransportVelocity")),
      limiter_(h_ref_ * h_ref_)
{
    static_assert(std::is_base_of<Limiter, LimiterType>::value,
                  "Limiter is not the base of LimiterType!");
    this->particles_->template addVariableToWrite<Vecd>("ExtraTransportVelocity");
}
//=================================================================================================//
template <class LimiterType, typename... CommonControlTypes>
void ExtraTransportForce<Inner<LimiterType>, CommonControlTypes...>::initialization(size_t index_i, Real dt)
{
    //** This initialisation is used to calculate the modified quantity of pos of the Transport Vel Correction(TVC) *
    //** the calculation of dr_tilde should be same as that in TVC */
    //** This step is redundant if we could modify the library to directly store and get dr_tilde */
    extra_transport_stress_[index_i] = Matd::Zero();
    if (this->within_scope_(index_i))
    {
        Real correction_scaling = 0.2 * h_ref_ * h_ref_;
        Real inv_h_ratio = 1.0;
        Real squared_norm = this->zero_gradient_residue_[index_i].squaredNorm();

        Vecd dr_tilde = correction_scaling * limiter_(squared_norm) * this->zero_gradient_residue_[index_i] * inv_h_ratio * inv_h_ratio;

        extra_transport_stress_[index_i] = this->rho_[index_i] * this->vel_[index_i] * dr_tilde.transpose();
    }
}
//=================================================================================================//
template <class LimiterType, typename... CommonControlTypes>
void ExtraTransportForce<Inner<LimiterType>, CommonControlTypes...>::interaction(size_t index_i, Real dt)
{
    extra_transport_vel_[index_i] = Vecd::Zero();
    if (this->within_scope_(index_i))
    {
        Matd extra_transport_stress_i = extra_transport_stress_[index_i];
        Vecd extra_transport_vel = Vecd::Zero();
        const Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd nablaW_ijV_j = inner_neighborhood.dW_ij_[n] * this->Vol_[index_j] * inner_neighborhood.e_ij_[n];

            extra_transport_vel += (extra_transport_stress_i + extra_transport_stress_[index_j]) * nablaW_ijV_j;
        }
        extra_transport_vel_[index_i] = extra_transport_vel / this->rho_[index_i];
    }
}
//=================================================================================================//
template <class LimiterType, typename... CommonControlTypes>
void ExtraTransportForce<Inner<LimiterType>, CommonControlTypes...>::update(size_t index_i, Real dt)
{
    if (this->within_scope_(index_i))
    {
        this->vel_[index_i] += extra_transport_vel_[index_i];
    }
}
//=================================================================================================//
template <typename... CommonControlTypes>
ExtraTransportForce<Contact<Boundary>, CommonControlTypes...>::ExtraTransportForce(BaseContactRelation &contact_relation)
    : ExtraTransportForce<Base, DataDelegateContact, CommonControlTypes...>(contact_relation),
      extra_transport_stress_(*this->particles_->template getVariableDataByName<Matd>("ExtraTransportStress")),
      extra_transport_vel_(*this->particles_->template getVariableDataByName<Vecd>("ExtraTransportVelocity"))
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        wall_Vol_.push_back(&(this->contact_particles_[k]->Vol_));
    }
}
//=================================================================================================//
template <typename... CommonControlTypes>
void ExtraTransportForce<Contact<Boundary>, CommonControlTypes...>::interaction(size_t index_i, Real dt)
{
    if (this->within_scope_(index_i))
    {
        Matd extra_transport_stress_i = extra_transport_stress_[index_i];
        Vecd extra_transport_vel = Vecd::Zero();
        for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
        {
            StdLargeVec<Real> &wall_Vol_k = *(wall_Vol_[k]);
            Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Vecd nablaW_ijV_j = contact_neighborhood.dW_ij_[n] * wall_Vol_k[index_j] * contact_neighborhood.e_ij_[n];

                extra_transport_vel += (extra_transport_stress_i)*nablaW_ijV_j;
            }
        }
        extra_transport_vel_[index_i] += extra_transport_vel / this->rho_[index_i];
    }
}

//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
TurbulentLinearGradientCorrectionMatrix<DataDelegationType>::
    TurbulentLinearGradientCorrectionMatrix(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      Vol_(*this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      turbu_B_(*this->particles_->template registerSharedVariable<Matd>(
          "TurbulentLinearGradientCorrectionMatrix", IdentityMatrix<Matd>::value)),
      B_(*this->particles_->template getVariableDataByName<Matd>("LinearGradientCorrectionMatrix"))
{
    this->particles_->template addVariableToWrite<Matd>("TurbulentLinearGradientCorrectionMatrix");
    this->particles_->template addVariableToSort<Matd>("TurbulentLinearGradientCorrectionMatrix");
    this->particles_->template addVariableToWrite<Matd>("LinearGradientCorrectionMatrix");
    this->particles_->template addVariableToSort<Matd>("LinearGradientCorrectionMatrix");
}

//=================================================================================================//

} // namespace fluid_dynamics
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//