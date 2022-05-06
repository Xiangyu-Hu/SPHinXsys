/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * --------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
 * and HU1527/12-1.															*
 *                                                                           *
 * Portions copyright (c) 2017-2020 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * --------------------------------------------------------------------------*/
/**
 * @file 	general_dynamics.h
 * @brief 	This is the particle dynamics aplliable for all type bodies
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef GENERAL_DYNAMICS_H
#define GENERAL_DYNAMICS_H

#include "all_particle_dynamics.h"
#include "base_body.h"
#include "base_particles.h"
#include "external_force.h"

#include <limits>

namespace SPH
{
	typedef DataDelegateSimple<SPHBody, BaseParticles, BaseMaterial> GeneralDataDelegateSimple;
	typedef DataDelegateInner<SPHBody, BaseParticles, BaseMaterial> GeneralDataDelegateInner;
	typedef DataDelegateContact<SPHBody, BaseParticles, BaseMaterial,
								SPHBody, BaseParticles, BaseMaterial, DataDelegateEmptyBase>
		GeneralDataDelegateContact;
	/**
	 * @class TimeStepInitialization
	 * @brief initialize a time step for a body.
	 * including initialize particle acceleration
	 * induced by viscous, gravity and other forces,
	 * set the number of ghost particles into zero.
	 */
	class TimeStepInitialization
		: public ParticleDynamicsSimple,
		  public GeneralDataDelegateSimple
	{
	private:
		UniquePtrKeeper<Gravity> gravity_ptr_keeper_;

	public:
		explicit TimeStepInitialization(SPHBody &sph_body);
		TimeStepInitialization(SPHBody &sph_body, Gravity &gravity);
		virtual ~TimeStepInitialization(){};

	protected:
		StdLargeVec<Vecd> &pos_n_, &dvel_dt_prior_;
		Gravity *gravity_;
		virtual void setupDynamics(Real dt = 0.0) override;
		virtual void Update(size_t index_i, Real dt = 0.0) override;
	};

	/**
	 * @class RandomizePartilePosition
	 * @brief Randomize the initial particle position
	 */
	class RandomizePartilePosition
		: public ParticleDynamicsSimple,
		  public GeneralDataDelegateSimple
	{
	public:
		explicit RandomizePartilePosition(SPHBody &sph_body);
		virtual ~RandomizePartilePosition(){};

	protected:
		StdLargeVec<Vecd> &pos_n_;
		Real randomize_scale_;
		virtual void Update(size_t index_i, Real dt = 0.0) override;
	};

	/**
	 * @class ParticleSmoothing
	 * @brief computing smoothed variable field by averaging with neighbors
	 */
	template <typename VariableType>
	class ParticleSmoothing : public InteractionDynamicsWithUpdate, public GeneralDataDelegateInner
	{
	public:
		explicit ParticleSmoothing(BaseBodyRelationInner &inner_relation, const std::string &variable_name)
			: InteractionDynamicsWithUpdate(*inner_relation.sph_body_),
			  GeneralDataDelegateInner(inner_relation),
			  W0_(body_->sph_adaptation_->getKernel()->W0(Vecd(0))),
			  smoothed_(*particles_->getVariableByName<VariableType>(variable_name))
		{
			particles_->registerAVariable(temp_, variable_name + "_temp");
		}

		virtual ~ParticleSmoothing(){};

	protected:
		const Real W0_;
		StdLargeVec<VariableType> &smoothed_, temp_;

		virtual void Interaction(size_t index_i, Real dt = 0.0) override
		{
			Real weight = W0_;
			VariableType summation = W0_ * smoothed_[index_i];
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				summation += inner_neighborhood.W_ij_[n] * smoothed_[index_j];
				weight += inner_neighborhood.W_ij_[n];
			}
			temp_[index_i] = summation / (weight + TinyReal);
		};

		virtual void Update(size_t index_i, Real dt = 0.0) override
		{
			smoothed_[index_i] = temp_[index_i];
		};
	};

	/**
	 * @class VelocityBoundCheck
	 * @brief  check whether particle velocity within a given bound
	 */
	class VelocityBoundCheck : public ParticleDynamicsReduce<bool, ReduceOR>,
							   public GeneralDataDelegateSimple
	{
	public:
		VelocityBoundCheck(SPHBody &sph_body, Real velocity_bound);
		virtual ~VelocityBoundCheck(){};

	protected:
		StdLargeVec<Vecd> &vel_n_;
		Real velocity_bound_;
		bool ReduceFunction(size_t index_i, Real dt = 0.0) override;
	};

	/**
	 * @class UpperFrontInXDirection
	 * @brief Get the upper front In X Direction for a SPH body
	 */
	class UpperFrontInXDirection : public ParticleDynamicsReduce<Real, ReduceMax>,
								   public GeneralDataDelegateSimple
	{
	public:
		explicit UpperFrontInXDirection(SPHBody &sph_body);
		virtual ~UpperFrontInXDirection(){};

	protected:
		StdLargeVec<Vecd> &pos_n_;
		Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
	};

	/**
	 * @class MaximumSpeed
	 * @brief Get the maximum particle speed in a SPH body
	 */
	class MaximumSpeed : public ParticleDynamicsReduce<Real, ReduceMax>,
						 public GeneralDataDelegateSimple
	{
	public:
		explicit MaximumSpeed(SPHBody &sph_body);
		virtual ~MaximumSpeed(){};

	protected:
		StdLargeVec<Vecd> &vel_n_;
		Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
	};

	/**
	 * @class BodyLowerBound
	 * @brief the lower bound of a body by reduced particle positions.
	 */
	class BodyLowerBound : public ParticleDynamicsReduce<Vecd, ReduceLowerBound>,
						   public GeneralDataDelegateSimple
	{
	public:
		explicit BodyLowerBound(SPHBody &sph_body);
		virtual ~BodyLowerBound(){};

	protected:
		StdLargeVec<Vecd> &pos_n_;
		Vecd ReduceFunction(size_t index_i, Real dt = 0.0) override;
	};

	/**
	 * @class BodyUpperBound
	 * @brief the upper bound of a body by reduced particle positions.
	 */
	class BodyUpperBound : public ParticleDynamicsReduce<Vecd, ReduceUpperBound>,
						   public GeneralDataDelegateSimple
	{
	public:
		explicit BodyUpperBound(SPHBody &sph_body);
		virtual ~BodyUpperBound(){};

	protected:
		StdLargeVec<Vecd> &pos_n_;
		Vecd ReduceFunction(size_t index_i, Real dt = 0.0) override;
	};

	/**
	 * @class BodySummation
	 * @brief Compute the summation of  a particle variable in a body
	 */
	template <typename VariableType>
	class BodySummation : public ParticleDynamicsReduce<VariableType, ReduceSum<VariableType>>,
						  public GeneralDataDelegateSimple
	{
	public:
		explicit BodySummation(SPHBody &sph_body, const std::string &variable_name)
			: ParticleDynamicsReduce<VariableType, ReduceSum<VariableType>>(sph_body),
			  GeneralDataDelegateSimple(sph_body),
			  variable_(*particles_->getVariableByName<VariableType>(variable_name))
		{
			this->initial_reference_ = VariableType(0);
		};
		virtual ~BodySummation(){};

	protected:
		StdLargeVec<VariableType> &variable_;
		VariableType ReduceFunction(size_t index_i, Real dt = 0.0) override
		{
			return variable_[index_i];
		};
	};

	/**
	 * @class BodyMoment
	 * @brief Compute the moment of a body
	 */
	template <typename VariableType>
	class BodyMoment : public BodySummation<VariableType>
	{
	public:
		explicit BodyMoment(SPHBody &sph_body, const std::string &variable_name)
			: BodySummation<VariableType>(sph_body, variable_name),
			  mass_(this->particles_->mass_){};
		virtual ~BodyMoment(){};

	protected:
		StdLargeVec<Real> &mass_;
		VariableType ReduceFunction(size_t index_i, Real dt = 0.0) override
		{
			return mass_[index_i] * this->variable_[index_i];
		};
	};

	/**
	 * @class TotalMechanicalEnergy
	 * @brief Compute the total mechanical (kinematic and potential) energy
	 */
	class TotalMechanicalEnergy
		: public ParticleDynamicsReduce<Real, ReduceSum<Real>>,
		  public GeneralDataDelegateSimple
	{
	private:
		UniquePtrKeeper<Gravity> gravity_ptr_keeper_;

	public:
		explicit TotalMechanicalEnergy(SPHBody &sph_body);
		TotalMechanicalEnergy(SPHBody &sph_body, Gravity &gravity_ptr);
		virtual ~TotalMechanicalEnergy(){};

	protected:
		StdLargeVec<Real> &mass_;
		StdLargeVec<Vecd> &vel_n_, &pos_n_;
		Gravity *gravity_;
		Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
	};
}
#endif // GENERAL_DYNAMICS_H