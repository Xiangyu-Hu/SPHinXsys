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
* @file    diffusion_splitting_base.h
* @brief   This is the splitting based method for diffusion optimization.
* @author  Bo Zhang and Xiangyu Hu
*/

#ifndef DIFFUSION_SPLITTING_BASE_H
#define DIFFUSION_SPLITTING_BASE_H

#include "all_particle_dynamics.h"
#include "diffusion_reaction.h"
#include "particle_dynamics_dissipation.h"

namespace SPH
{
	/**
	 * @class OptimizationBySplittingAlgorithmBase
	 * @brief The base class for optimization using the splitting algorithm.
	 */
	template <class ParticlesType, typename VariableType>
	class OptimizationBySplittingAlgorithmBase
		: public LocalDynamics, public DataDelegateInner<ParticlesType>
	{
	public:
		explicit OptimizationBySplittingAlgorithmBase(BaseInnerRelation& inner_relation, 
			                                          const std::string &variable_name);
		virtual ~OptimizationBySplittingAlgorithmBase() {};

	public:
		size_t phi_;
		StdLargeVec<Real>& Vol_, & mass_;
		StdLargeVec<Vecd>& normal_vector_;
		StdLargeVec<VariableType>& variable_;
		StdLargeVec<Real>& heat_flux_, & heat_source_;
		StdLargeVec<int> &splitting_index_;
		StdLargeVec<Real> &species_modified_, &species_recovery_;
		StdLargeVec<Real> &parameter_recovery_, &eta_regularization_;
		StdLargeVec<Real> &residual_T_local_, &residual_T_global_;
		StdLargeVec<Real> &residual_k_local_, &residual_k_global_;
		StdLargeVec<Real> &variation_local_, &variation_global_;
		StdLargeVec<Real> &residual_after_splitting_;
       
		
		StdVec<BaseDiffusion*> all_diffusion_;

		virtual void interaction(size_t index_i, Real dt = 0.0) = 0;
	};

	/**
	 * @class RegularizationByDiffusion
	 * @brief Regularize the optimized parameter by diffusion analogy method 
	 *        after each splitting step, which could smooth the distribution 
	 *        and avoid some local optimal solution.
	 */
	template <class ParticlesType, typename VariableType>
	class RegularizationByDiffusionAnalogy
		: public OptimizationBySplittingAlgorithmBase<ParticlesType, VariableType>
    {
	public:
        RegularizationByDiffusionAnalogy(BaseInnerRelation &inner_relation, const std::string &variable_name,
                                         Real initial_eta = 1, Real variation = 1);
        virtual ~RegularizationByDiffusionAnalogy(){};

		void UpdateCurrentEta(Real initial_eta) { initial_eta_ = initial_eta; }
		void UpdateMaximumVariation(Real maximum_variation) { maximum_variation_ = maximum_variation; }
		void UpdateAverageVariation(Real averaged_variation) { averaged_variation_ = averaged_variation; }

	protected:
		Real initial_eta_, maximum_variation_, averaged_variation_;
		virtual ErrorAndParameters<VariableType> computeVariationAndParameters(size_t index_i, Real dt);
		virtual void updateStatesByVariation(size_t index_i, Real dt, const ErrorAndParameters<VariableType>& variation_and_parameters);
		virtual void interaction(size_t index_i, Real dt = 0.0) override;
    };

	/**
	 * @class UpdateRegularizationVariation
	 * @brief The global variation of parameter can be updated after the process 
	 *        of splitting which is an essential parameter for optimization schedule.
	 */
	template <class ParticlesType, typename VariableType>
    class UpdateRegularizationVariation
        : public OptimizationBySplittingAlgorithmBase<ParticlesType, VariableType>
    {
      public:
        UpdateRegularizationVariation(BaseInnerRelation &inner_relation, const std::string &variable_name);
        virtual ~UpdateRegularizationVariation(){};

      protected:
        /* Redefine the compute function to avoid non-meaningful initial variation. */
        virtual ErrorAndParameters<VariableType> computeVariationAndParameters(size_t index_i, Real dt);
        virtual void interaction(size_t index_i, Real dt = 0.0) override;
    };
}
#endif //DIFFUSION_SPLITTING_BASE_H