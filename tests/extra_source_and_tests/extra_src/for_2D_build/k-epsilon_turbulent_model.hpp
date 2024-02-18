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
	BaseTurtbulentModel<Base, DataDelegationType>::BaseTurtbulentModel(BaseRelationType& base_relation)
		: LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
		particle_spacing_min_(base_relation.getSPHBody().sph_adaptation_->MinimumSpacing()),
		rho_(this->particles_->rho_), vel_(this->particles_->vel_),
		mu_(DynamicCast<Fluid>(this, this->particles_->getBaseMaterial()).ReferenceViscosity()), dimension_(Vecd(0).size()),
		smoothing_length_(this->sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
		turbu_k_(*(this->particles_->template registerSharedVariable<Real>("TurbulenceKineticEnergy"))),
		turbu_k_prior_(*(this->particles_->template registerSharedVariable<Real>("TurbulenceKineticEnergyPrior"))),
		turbu_epsilon_(*(this->particles_->template registerSharedVariable<Real>("TurbulentDissipation"))),
		turbu_epsilon_prior_(*(this->particles_->template registerSharedVariable<Real>("TurbulentDissipationPrior"))),
		turbu_mu_(*(this->particles_->template registerSharedVariable<Real>("TurbulentViscosity")))
	{}
//=================================================================================================//
	template <class DataDelegationType>
	template <class BaseRelationType>
	TKEnergyForce<Base, DataDelegationType>::
		TKEnergyForce(BaseRelationType& base_relation) :
		BaseTurtbulentModel<Base, DataDelegationType>(base_relation),
		force_(this->particles_->force_), mass_(this->particles_->mass_),
		indicator_(this->particles_->indicator_), pos_(this->particles_->pos_),
		turbu_k_(*this->particles_->template getVariableByName<Real>("TurbulenceKineticEnergy")),
		test_k_grad_rslt_(*this->particles_->template registerSharedVariable<Vecd>("TkeGradResult")){}
//=================================================================================================//
	template <class DataDelegationType>
	template <class BaseRelationType>
	GetVelocityGradient<Base, DataDelegationType>::
		GetVelocityGradient(BaseRelationType& base_relation) 
		:LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
		vel_(this->particles_->vel_), pos_(this->particles_->pos_),
		is_near_wall_P1_(*this->particles_->template getVariableByName<int>("IsNearWallP1"))
	{
		this->particles_->registerVariable(velocity_gradient_, "VelocityGradient");
		this->particles_->registerSortableVariable<Matd>("VelocityGradient");
		this->particles_->addVariableToWrite<Matd>("VelocityGradient");

		//for test
		this->particles_->registerVariable(velocity_gradient_wall, "Velocity_Gradient_Wall");
		this->particles_->registerSortableVariable<Matd>("Velocity_Gradient_Wall");
		this->particles_->addVariableToWrite<Matd>("Velocity_Gradient_Wall");
	}
//=================================================================================================//
	template <class DataDelegationType>
	template <class BaseRelationType>
	TurbuViscousForce<DataDelegationType>::TurbuViscousForce(BaseRelationType& base_relation)
		: ViscousForce<DataDelegationType>(base_relation),
		turbu_k_(*this->particles_->template getVariableByName<Real>("TurbulenceKineticEnergy")),
		turbu_mu_(*this->particles_->template getVariableByName<Real>("TurbulentViscosity")),
		wall_Y_plus_(*this->particles_->template getVariableByName<Real>("WallYplus")),
		wall_Y_star_(*this->particles_->template getVariableByName<Real>("WallYstar")),
		velo_friction_(*this->particles_->template getVariableByName<Vecd>("FrictionVelocity")),
		y_p_(*this->particles_->template getVariableByName<Real>("Y_P")),
		molecular_viscosity_(DynamicCast<Fluid>(this, this->particles_->getBaseMaterial()).ReferenceViscosity()) {}
//=================================================================================================//

}
//=================================================================================================//
}
//=================================================================================================//