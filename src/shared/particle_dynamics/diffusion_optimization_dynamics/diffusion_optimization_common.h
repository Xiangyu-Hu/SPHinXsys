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
* @file 	diffusion_optimization_common.h
* @brief 	This is the method for updating, calculating error, observed parameters.
* @author   Bo Zhang and Xiangyu Hu
*/

#ifndef DIFFUSION_OPTIMIZATION_COMMON_H
#define DIFFUSION_OPTIMIZATION_COMMON_H

#include "diffusion_splitting_base.h"
#include "diffusion_splitting_parameter.h"
#include "diffusion_splitting_state.h"
#include "all_particle_dynamics.h"

namespace SPH
{
	/**
	 * @class UpdateUnitNormalVector
	 * @brief Initialize and update the unit normal vector to boundary.
	 */
	template <class BaseParticlesType, class BaseMaterialType, class ContactBaseParticlesType, 
		      class ContactBaseMaterialType, int NUM_SPECIES = 1>
	class UpdateUnitNormalVector
		: public LocalDynamics,
		  public DiffusionReactionInnerData<BaseParticlesType, BaseMaterialType>,
		  public DiffusionReactionContactData<BaseParticlesType, BaseMaterialType, 
		                                      ContactBaseParticlesType, ContactBaseMaterialType>
	{
	protected:
		StdLargeVec<Vecd>& normal_vector_;

	public:
		UpdateUnitNormalVector(ComplexRelation& body_complex_relation);
		virtual ~UpdateUnitNormalVector() {};
		virtual void interaction(size_t index_i, Real dt = 0.0);
	};

	/**
	 * @class UpdateNormalDistance
	 * @brief Calculate the distance to the boundary.
	 */
	template<class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES = 1>
	class UpdateNormalDistance
		: public LocalDynamics,
		  public DiffusionReactionInnerData< BaseParticlesType, BaseMaterialType, NUM_SPECIES>
	{
	protected:
		Real W0_, sigma0_, cutoff_radius_;
		StdLargeVec<Real>& normal_distance_;

	public:
		UpdateNormalDistance(BaseInnerRelation& inner_relation);
		virtual ~UpdateNormalDistance() {};
		virtual void interaction(size_t index_i, Real dt = 0.0);
	};

	/**
	 * @class  ComputeAveragedErrorOrPositiveParameter
	 * @brief  Computing the averaged error on the (partly) optimization domain
	           or any parameter only with positive values.
	 */
	template <class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES = 1>
	class ComputeAveragedErrorOrPositiveParameter
		: public DiffusionReactionSpeciesSummation<BaseParticlesType, BaseMaterialType, NUM_SPECIES>
	{
	protected:
		StdLargeVec<Real>& variable_;
		
	public:
		ComputeAveragedErrorOrPositiveParameter(SPHBody &sph_body, const std::string &error_name)
			: DiffusionReactionSpeciesSummation<BaseParticlesType, BaseMaterialType,
			                                    NUM_SPECIES>(sph_body, error_name),
			  variable_(*this->particles_->template getVariableByName<Real>(error_name)) {};
		ComputeAveragedErrorOrPositiveParameter(BodyPartByParticle& body_part, const std::string &error_name)
			: ComputeAveragedErrorOrPositiveParameter(body_part.getSPHBody(), error_name) {};
		virtual ~ComputeAveragedErrorOrPositiveParameter() {};

		Real reduce(size_t index_i, Real dt =0.0)
		{
			return abs(variable_[index_i]);
		}	  
	};

	/*
	 * @Class ComputeMaximumError
	 * @brief get and return the maximum PDE residual among the whole domain.
	 */
	template <class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES = 1>
	class ComputeMaximumError 
		: public LocalDynamicsReduce<Real, ReduceMax>,
		  public DiffusionReactionSimpleData<BaseParticlesType, BaseMaterialType, NUM_SPECIES>
	{
	public:
		ComputeMaximumError(SPHBody& sph_body, const std::string& error_name)
			: LocalDynamicsReduce<Real, ReduceMax>(sph_body, Real(0)),
			  DiffusionReactionSimpleData<BaseParticlesType, BaseMaterialType, NUM_SPECIES>(sph_body),
			  variable_(*this->particles_->template getVariableByName<Real>(error_name)) {};
		ComputeMaximumError(BodyPartByParticle& body_part, const std::string& error_name)
			: ComputeMaximumError(body_part.getSPHBody(), error_name) {};
		virtual ~ComputeMaximumError() {};

	protected:
		StdLargeVec<Real>& variable_;
		Real reduce(size_t index_i, Real dt = 0.0)
		{
			return abs(variable_[index_i]);
		};
	};

	/**
	 * @class ThermalDiffusivityConstrain
	 * @brief The thermal diffusivity on each particle will be corrected with
	 *		  the same ratio according to the total thermal diffusivity.
	 */
	template<class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES = 1>
	class ThermalDiffusivityConstrain
		: public LocalDynamics,
		  public DiffusionReactionSimpleData<BaseParticlesType, BaseMaterialType, NUM_SPECIES>
	{
	public:
		ThermalDiffusivityConstrain(SPHBody& diffusion_body, const std::string& variable_name,
			                        Real initial_thermal_diffusivity = 1);
		virtual ~ThermalDiffusivityConstrain() {};
		void UpdateAveragedParameter(Real new_averaged_thermal_diffusivity);
	
	protected:
		Real initial_thermal_diffusivity_;
		Real new_averaged_thermal_diffusivity_;
		StdLargeVec<VariableType>& local_thermal_conductivity_;
		void update(size_t index_i, Real dt = 0.0);
	};
}

#endif DIFFUSION_OPTIMIZATION_COMMON_H