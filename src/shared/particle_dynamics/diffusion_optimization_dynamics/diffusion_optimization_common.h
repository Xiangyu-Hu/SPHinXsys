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
	 * @class  ComputeAverageErrorOrPositiveParameter
	 * @breif  Computing the average error on the (partly) optimization 
	 *         domain or any paraeter only with positive values.
	 */
	template <class DynamicsIdentifier, class ParticlesType>
	class ComputeAverageErrorOrPositiveParameter
		: public SpeciesSummation<DynamicsIdentifier, ParticlesType>
	{
    protected:
		StdLargeVec<Real> &variable_;

	public:
		ComputeAverageErrorOrPositiveParameter(DynamicsIdentifier &identifier, const std::string &variable_name):
			SpeciesSummation<DynamicsIdentifier, ParticlesType>(identifier, variable_name),
			variable_(*this->particles_->template getVariableByName<Real>(variable_name)){};
		virtual ~ComputeAverageErrorOrPositiveParameter(){};
		
		Real reduce(size_t index_i, Real dt = 0.0)
		{
			return abs(variable_[index_i]);
		};
	};

	/**
	 * @class ComputeMaximumError
	 * @breif Get and return the maximum residual among the whole domain
	 */
	template <class DynamicsIdentifier, class ParticlesType>
    class ComputeMaximumError
		: public BaseLocalDynamicsReduce<Real, ReduceMax, DynamicsIdentifier>,
          public DiffusionReactionSimpleData<ParticlesType>
    {
    protected:
		StdLargeVec<Real> &variable_;

	public:
		ComputeMaximumError(DynamicsIdentifier &identifier, const std::string &variable_name)
			: BaseLocalDynamicsReduce<Real, ReduceMax, DynamicsIdentifier>(identifier, Real(0)),
              DiffusionReactionSimpleData<ParticlesType>(identifier.getSPHBody()),
              variable_(*this->particles_->template getVariableByName<Real>(variable_name)){};

		Real reduce(size_t index_i, Real dt = 0.0)
        {
			return abs(variable_[index_i]);
        }
    };

	/**
	 * @class ThermalConductivityConstrain
	 * @brief The thermal diffusivity on each particle will be corrected with 
	 *        the same ratio according to the total thermal diffusivity.
	 */
	template <class ParticlesType>
    class ThermalConductivityConstrain
        : public LocalDynamics,
          public DiffusionReactionSimpleData<ParticlesType>
    {
    public:
        ThermalConductivityConstrain(SPHBody &diffusion_body, const std::string &variable_name,
                                     Real initial_thermal_conductivity = 1);
        virtual ~ThermalConductivityConstrain(){};
        void UpdateAverageParameter(Real new_average_thermal_diffusivity);

    protected:
        Real initial_thermal_conductivity_;
        Real new_average_thermal_conductivity_;
        StdLargeVec<Real> &local_thermal_conductivity_;
        void update(size_t index_i, Real dt = 0.0);
    };
}

#endif //DIFFUSION_OPTIMIZATION_COMMON_H