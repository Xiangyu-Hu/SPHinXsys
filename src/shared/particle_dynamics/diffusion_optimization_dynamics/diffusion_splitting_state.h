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
* @file 	diffusion_splitting_state.h
* @brief 	This is the splitting method for solving state field in optimization problem.
* @author   Bo Zhang and Xiangyu H
*/

#ifndef DIFFUSION_SPLITTING_STATE_H
#define DIFFUSION_SPLITTING_STATE_H

#include "diffusion_splitting_base.h"
#include "diffusion_splitting_parameter.h"

namespace SPH
{
	/**
	 * @class TemperatureSplittingByPDEInner
	 * @brief The temperature on each particle will be modified innerly to satisfy the PDEs.
	 */
	template <class ParticlesType, typename VariableType>
	class TemperatureSplittingByPDEInner
		: public OptimizationBySplittingAlgorithmBase<ParticlesType, VariableType>
	{
	public:
		TemperatureSplittingByPDEInner(BaseInnerRelation &inner_relation, const std::string &variable_name);
		virtual ~TemperatureSplittingByPDEInner(){};
	
	protected:
		virtual ErrorAndParameters<VariableType> computeErrorAndParameters(size_t index_i, Real dt = 0.0);
		virtual void updateStatesByError(size_t index_i, Real dt, const ErrorAndParameters<VariableType> &error_and_parameters);
		virtual void interaction(size_t index_i, Real dt = 0.0) override;
};

	/**
	 * @class TemperatureSplittingByPDEWithBoundary
	 * @brief The temperature on each particle will be modified with boundary to satisfy the PDEs.
	 */
	template <class ParticlesType, class ContactParticlesType, typename VariableType>
	class TemperatureSplittingByPDEWithBoundary
		: public TemperatureSplittingByPDEInner<ParticlesType, VariableType>,
		  public DataDelegateContact<ParticlesType, ContactParticlesType, DataDelegateEmptyBase>
	{
	public:
		TemperatureSplittingByPDEWithBoundary(ComplexRelation &complex_relation, const std::string &variable_name);
		virtual ~TemperatureSplittingByPDEWithBoundary(){};

	protected:
		StdVec<StdLargeVec<VariableType> *> boundary_variable_;
		StdVec<StdLargeVec<Real> *> boundary_heat_flux_;
		StdVec<StdLargeVec<Vecd> *> boundary_normal_vector_;
		virtual ErrorAndParameters<VariableType> computeErrorAndParameters(size_t index_i, Real dt = 0.0) override;
	};

	/**
	 * @class UpdateTemperaturePDEResidual
	 * @brief Update the global residual from temperature after one splitting loop finished.
	 */
	template <typename TemperatureSplittingType, typename BaseBodyRelationType, typename VariableType>
        class UpdateTemperaturePDEResidual : public TemperatureSplittingType
    {
	public:
		UpdateTemperaturePDEResidual(BaseBodyRelationType& body_relation, const std::string& variable_name);
		virtual ~UpdateTemperaturePDEResidual() {};

	protected:
		virtual void interaction(size_t index_i, Real dt = 0.0) override;
    };
}

#endif //DIFFUSION_SPLITTING_STATE_H