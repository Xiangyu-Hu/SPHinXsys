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
* @file 	diffusion_splitting_parameter.h
* @brief    This is the method of splitting parameter for reversing or optimization problem.
* @author   Bo Zhang and Xiangyu Hu
*/

#ifndef DIFFUSION_SPLITTING_PARAMETER_H
#define DIFFUSION_SPLITTING_PARAMETER_H

#include "diffusion_splitting_base.h"
#include "diffusion_splitting_state.h"

namespace SPH
{
	/* Here the parameter splitting may use the pseudo time-step to keep with
	   temperature splitting, and the diagonal dominance will stable the calculation.*/

	/**
	  * @class ParameterSplittingByPDEInner
	  * @brief Modify the parameter inner by splitting operator based on PDE.
	  */
	template <class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES = 1>
	class ParameterSplittingByPDEInner
		: public OptimizationBySplittingAlgorithmBase<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>
	{
	public:
		ParameterSplittingByPDEInner(BaseInnerRelation& inner_relation, const std::string& variable_name);
		virtual ~ParameterSplittingByPDEInner() {};

	protected:
		virtual ErrorAndParameters<VariableType> computeErrorAndParameters(size_t index_i, Real dt = 0.0);
		virtual void updateStatesByError(size_t index_i, Real dt, const ErrorAndParameters<VariableType>& error_and_parameters);
		virtual void interaction(size_t index_i, Real dt = 0.0) override;
	};

	/**
	  * @class ParameterSplittingByPDEWithBoundary
	  * @brief Modify the parameter contact with the boundary by splitting operator based on PDE.
	  */
	template <class BaseParticlesType, class BaseMaterialType, class ContactBaseParticlesType, 
		      class ContactBaseMaterialType, typename VariableType, int NUM_SPECIES = 1>
	class ParameterSplittingByPDEWithBoundary
		: public ParameterSplittingByPDEInner<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>,
		  public DiffusionReactionContactData<BaseParticlesType, BaseMaterialType, ContactBaseParticlesType,
		                                      ContactBaseMaterialType, NUM_SPECIES>
	{
	public:
		ParameterSplittingByPDEWithBoundary(ComplexRelation& complex_relation, const std::string& variable_name);
		virtual ~ParameterSplittingByPDEWithBoundary() {};

	protected:
		StdVec<StdLargeVec<Vecd>*> boundary_normal_vector_;
		StdVec<StdLargeVec<Real>*> boundary_heat_flux_, boundary_normal_distance_;
		StdVec<StdVec<StdLargeVec<Real>>*> boundary_species_;
		virtual ErrorAndParameters<VariableType> computeErrorAndParameters(size_t index_i, Real dt = 0.0) override;
	};
		                   
	/**
	  * @class UpdateParameterPDEResidual
	  * @brief Update the global residual from parameter after one splitting loop finished.
	  */
    template <typename ParameterSplittingType, typename BaseBodyRelationType, typename VariableType>
	class UpdateParameterPDEResidual : public ParameterSplittingType
	{
	public:
		UpdateParameterPDEResidual(BaseBodyRelationType &body_relation, const std::string &variable_name);
		virtual ~UpdateParameterPDEResidual() {};

	protected:
		virtual void interaction(size_t index_i, Real dt = 0.0) override;
	};

	/**
	  * @class ParameterSplittingByBCInner
	  * @brief Modify the parameter inner by splitting operator based on BC constraint.
	  *        Used for the boundary particle is within the thermal domain.
	  */
	/*template <class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES = 1>
	class ParameterSplittingByBCInner
		: public OptimizationBySplittingAlgorithmBase<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>
	{
	public:
		ParameterSplittingByBCInner(BaseInnerRelation& inner_relation, const std::string& variable_name);
		virtual ~ParameterSplittingByBCInner() {};

	protected:
		virtual ErrorAndParameters<VariableType> computeErrorAndParameters(size_t index_i, Real dt = 0.0);
		virtual void updateStatesByError(size_t index_i, Real dt, const ErrorAndParameters<VariableType>& error_and_parameters);
		virtual void interaction(size_t index_i, Real dt = 0.0) override;
	};*/

	/**
	  * @class ParameterSplittingByBCWithBoundary
	  * @brief Modify the parameter contact by splitting operator based on BC constraint.
	  *        Used for the boundary particle is in another body, i.e., dummy particle.
	  */
	/*template <class BaseParticlesType, class BaseMaterialType, class ContactBaseParticlesType,
		      class ContactBaseMaterialType, typename VariableType, int NUM_SPECIES = 1>
	class ParameterSplittingByBCWithBoundary 
		: public ParameterSplittingByBCInner<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>,
		  public DiffusionReactionContactData<BaseParticlesType, BaseMaterialType, 
		                                      ContactBaseParticlesType, ContactBaseMaterialType, NUM_SPECIES>
	{
	public:
		ParameterSplittingByBCWithBoundary(ComplexRelation& complex_relation, const std::string& variable_name);
		virtual ~ParameterSplittingByBCWithBoundary() {};

	protected:
		StdVec<StdLargeVec<int>*> boundary_boundary_index_;
		StdVec<StdLargeVec<Vecd>*> boundary_normal_vector_;
		StdVec<StdLargeVec<Real>*> boundary_heat_flux_, boundary_normal_distance_, boundary_parameter_recovery_;
		StdVec<StdLargeVec<VariableType>*> boundary_variable_;
		StdVec<StdVec<StdLargeVec<Real>>*> boundary_species_;
		virtual ErrorAndParameters<VariableType> computeErrorAndParameters(size_t index_i, Real dt = 0.0);
		virtual void updateStatesByError(size_t index_i, Real dt, const ErrorAndParameters<VariableType>& error_and_parameters);
		virtual void interaction(size_t index_i, Real dt = 0.0) override;
	};*/
	
	/**
	 * @class UpdateParameterBCResidual
	 * @brief Update the global residual from parameter after one splitting loop finished.
	 */
	/*template <typename ParameterSplittingType, typename BaseBodyRelationType, typename VariableType>
	class UpdateParameterBCResidual : public ParameterSplittingType
	{
	public:
		UpdateParameterBCResidual(BaseBodyRelationType& body_relation, const std::string& variable_name);
		virtual ~UpdateParameterBCResidual() {};

	protected:
		virtual void interaction(size_t index_i, Real dt = 0.0) override;
	};*/
}

#endif DIFFUSION_SPLITTING_PARAMETER_H

