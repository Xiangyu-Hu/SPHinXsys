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
BaseTurbulentModel<Base, DataDelegationType>::BaseTurbulentModel(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      turbu_strain_rate_(this->particles_->template registerStateVariable<Matd>("TurbulentStrainRate")),
      viscosity_(DynamicCast<Viscosity>(this, this->particles_->getBaseMaterial())),
      mu_(viscosity_.ReferenceViscosity()),
      smoothing_length_(this->sph_body_.getSPHAdaptation().ReferenceSmoothingLength()),
      particle_spacing_min_(base_relation.getSPHBody().getSPHAdaptation().MinimumSpacing()),
      rho_(this->particles_->template getVariableDataByName<Real>("Density")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      vel_(this->particles_->template getVariableDataByName<Vecd>("Velocity")),
      dimension_(2) {}

//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
TKEnergyForce<Base, DataDelegationType>::
    TKEnergyForce(BaseRelationType &base_relation)
    : BaseTurbulentModel<Base, DataDelegationType>(base_relation),
      force_(this->particles_->template registerStateVariable<Vecd>("Force")),
      mass_(this->particles_->template getVariableDataByName<Real>("Mass")),
      indicator_(this->particles_->template getVariableDataByName<int>("Indicator")),
      pos_(this->particles_->template getVariableDataByName<Vecd>("Position")),
      turbu_k_(this->particles_->template getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      test_k_grad_rslt_(this->particles_->template registerStateVariable<Vecd>("TkeGradResult")) {}
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
GetVelocityGradient<DataDelegationType>::
    GetVelocityGradient(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      vel_(this->particles_->template getVariableDataByName<Vecd>("Velocity")),
      pos_(this->particles_->template getVariableDataByName<Vecd>("Position")),
      is_near_wall_P1_(this->particles_->template getVariableDataByName<int>("IsNearWallP1")),
      is_near_wall_P2_(this->particles_->template getVariableDataByName<int>("IsNearWallP2")),
      velocity_gradient_(this->particles_->template registerStateVariable<Matd>("TurbulentVelocityGradient")),
      velocity_gradient_wall(this->particles_->template registerStateVariable<Matd>("Velocity_Gradient_Wall")) {}
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
TurbuViscousForce<DataDelegationType>::TurbuViscousForce(BaseRelationType &base_relation)
    : ViscousForce<DataDelegationType>(base_relation),
      turbu_k_(this->particles_->template getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_mu_(this->particles_->template getVariableDataByName<Real>("TurbulentViscosity")),
      wall_Y_plus_(this->particles_->template getVariableDataByName<Real>("WallYplus")),
      wall_Y_star_(this->particles_->template getVariableDataByName<Real>("WallYstar")),
      velo_friction_(this->particles_->template getVariableDataByName<Vecd>("FrictionVelocity")),
      y_p_(this->particles_->template getVariableDataByName<Real>("Y_P")),
      is_near_wall_P2_(this->particles_->template getVariableDataByName<int>("IsNearWallP2")),
      viscosity_(DynamicCast<Viscosity>(this, this->particles_->getBaseMaterial())),
      molecular_viscosity_(viscosity_.ReferenceViscosity()),
      c0_(DynamicCast<Fluid>(this, this->particles_->getBaseMaterial()).ReferenceSoundSpeed()) {}
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
TurbulentLinearGradientCorrectionMatrix<DataDelegationType>::
    TurbulentLinearGradientCorrectionMatrix(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      turbu_B_(this->particles_->template registerStateVariable<Matd>(
          "TurbulentLinearGradientCorrectionMatrix", IdentityMatrix<Matd>::value)),
      B_(this->particles_->template getVariableDataByName<Matd>("LinearGradientCorrectionMatrix"))
{
    this->particles_->template addVariableToSort<Matd>("TurbulentLinearGradientCorrectionMatrix");
    this->particles_->template addVariableToSort<Matd>("LinearGradientCorrectionMatrix");
}

//=================================================================================================//
} // namespace fluid_dynamics
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//