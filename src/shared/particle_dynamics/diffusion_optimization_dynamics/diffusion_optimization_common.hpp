/**
 * @file 	diffusion_optimization_common.hpp
 * @author	Bo Zhang and Xiangyu Hu
 */

#ifndef DIFFUSION_OPTIMIZATION_COMMON_HPP
#define DIFFUSION_OPTIMIZATION_COMMON_HPP

#include "diffusion_optimization_common.h"

namespace SPH
{
//=================================================================================================//
template <class DynamicsIdentifier>
ThermalConductivityConstraint<DynamicsIdentifier>::
    ThermalConductivityConstraint(DynamicsIdentifier &identifier, const std::string &variable_name,
                                  Real initial_thermal_conductivity)
    : LocalDynamics(identifier.getSPHBody()),

      initial_thermal_conductivity_(initial_thermal_conductivity),
      new_average_thermal_conductivity_(0.0),
      local_diffusivity_(this->particles_->template getVariableDataByName<Real>(variable_name)){};
//=================================================================================================//
template <class DynamicsIdentifier>
void ThermalConductivityConstraint<DynamicsIdentifier>::
    UpdateAverageParameter(Real new_average_thermal_conductivity)
{
    new_average_thermal_conductivity_ = new_average_thermal_conductivity;
};
//=================================================================================================//
template <class DynamicsIdentifier>
void ThermalConductivityConstraint<DynamicsIdentifier>::
    update(size_t index_i, Real dt)
{
    local_diffusivity_[index_i] = local_diffusivity_[index_i] *
                                  initial_thermal_conductivity_ / new_average_thermal_conductivity_;
}
//=================================================================================================//
} // namespace SPH
#endif // DIFFUSION_OPTIMIZATION_COMMON_HPP