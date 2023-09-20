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
	template <class ParticlesType>
	ThermalConductivityConstrain<ParticlesType>::
		ThermalConductivityConstrain(SPHBody &diffusion_body, const std::string &variable_name,
                                     Real initial_thermal_conductivity) : 
	    LocalDynamics(diffusion_body),                                    
		DiffusionReactionSimpleData<ParticlesType>(diffusion_body),                                                                  
		initial_thermal_conductivity_(initial_thermal_conductivity),                            
		new_average_thermal_conductivity_(0.0),                                                                 
		local_thermal_conductivity_(*this->particles_->template getVariableByName<Real>(variable_name)){};
	//=================================================================================================//
	template <class ParticlesType>
    void ThermalConductivityConstrain<ParticlesType>::
		UpdateAverageParameter(Real new_average_thermal_conductivity)
    {
		new_average_thermal_conductivity_ = new_average_thermal_conductivity;
    };
	//=================================================================================================//
	template <class ParticlesType>
    void ThermalConductivityConstrain<ParticlesType>::
		update(size_t index_i, Real dt)
	{
		local_thermal_conductivity_[index_i] = local_thermal_conductivity_[index_i] *
			initial_thermal_conductivity_ / new_average_thermal_conductivity_;
	}
	//=================================================================================================//
}
#endif //DIFFUSION_OPTIMIZATION_COMMON_HPP