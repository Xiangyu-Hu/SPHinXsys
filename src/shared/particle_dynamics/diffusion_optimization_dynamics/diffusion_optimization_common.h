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
 * @file diffusion_optimization_common.h
 * @brief This is the method for updating, calculating error, observed parameters.
 * TODO: The methods can be replaced by the general reduce methods in the future.
 * @author Bo Zhang and Xiangyu Hu
 */

#ifndef DIFFUSION_OPTIMIZATION_COMMON_H
#define DIFFUSION_OPTIMIZATION_COMMON_H

#include "all_particle_dynamics.h"
#include "diffusion_splitting_base.h"
#include "diffusion_splitting_parameter.h"
#include "diffusion_splitting_state.h"

namespace SPH
{
/**
 * @class  ComputeTotalErrorOrPositiveParameter
 * @brief  Computing the average error on the (partly) optimization
 *         domain or any parameter only with positive values.
 */
template <class DynamicsIdentifier>
class ComputeTotalErrorOrPositiveParameter
    : public BaseLocalDynamicsReduce<ReduceSum<Real>, DynamicsIdentifier>
{
  protected:
    Real *variable_;

  public:
    ComputeTotalErrorOrPositiveParameter(DynamicsIdentifier &identifier, const std::string &variable_name)
        : BaseLocalDynamicsReduce<ReduceSum<Real>, DynamicsIdentifier>(identifier),
          variable_(this->particles_->template getVariableDataByName<Real>(variable_name)){};
    virtual ~ComputeTotalErrorOrPositiveParameter(){};

    Real reduce(size_t index_i, Real dt = 0.0)
    {
        return abs(variable_[index_i]);
    };
};

/**
 * @class ComputeMaximumError
 * @brief Get and return the maximum residual among the whole domain
 */
template <class DynamicsIdentifier>
class ComputeMaximumError
    : public BaseLocalDynamicsReduce<ReduceMax, DynamicsIdentifier>
{
  protected:
    Real *variable_;

  public:
    ComputeMaximumError(DynamicsIdentifier &identifier, const std::string &variable_name)
        : BaseLocalDynamicsReduce<ReduceMax, DynamicsIdentifier>(identifier),
          variable_(this->particles_->template getVariableDataByName<Real>(variable_name)){};

    Real reduce(size_t index_i, Real dt = 0.0)
    {
        return abs(variable_[index_i]);
    }
};

/**
 * @class ThermalConductivityConstraint
 * @brief The thermal diffusivity on each particle will be corrected with
 *        the same ratio according to the total thermal diffusivity.
 */
template <class DynamicsIdentifier>
class ThermalConductivityConstraint
    : public LocalDynamics
{
  public:
    ThermalConductivityConstraint(DynamicsIdentifier &identifier, const std::string &variable_name,
                                  Real initial_thermal_conductivity = 1);
    virtual ~ThermalConductivityConstraint(){};
    void UpdateAverageParameter(Real new_average_thermal_diffusivity);

  protected:
    Real initial_thermal_conductivity_;
    Real new_average_thermal_conductivity_;
    Real *local_diffusivity_;
    void update(size_t index_i, Real dt = 0.0);
};
} // namespace SPH

#endif // DIFFUSION_OPTIMIZATION_COMMON_H