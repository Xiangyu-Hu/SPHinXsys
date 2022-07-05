/**
 * @file diffusion_reaction.cpp
 * @brief These are classes for diffusion and reaction properties
 * @author Chi Zhang and Xiangyu Hu
 */

#include "diffusion_reaction.h"

#include "base_particles.hpp"

namespace SPH
{
	//=================================================================================================//
	void DirectionalDiffusion::initializeDirectionalDiffusivity(Real diff_cf, Real bias_diff_cf, Vecd bias_direction)
	{
		bias_diff_cf_ = bias_diff_cf;
		bias_direction_ = bias_direction;
		Matd diff_i = diff_cf_ * Matd(1.0) + bias_diff_cf_ * SimTK::outer(bias_direction_, bias_direction_);
		transformed_diffusivity_ = inverseCholeskyDecomposition(diff_i);
	};
	//=================================================================================================//
	void LocalDirectionalDiffusion::assignBaseParticles(BaseParticles *base_particles)
	{
		DirectionalDiffusion::assignBaseParticles(base_particles);
		initializeFiberDirection();
	};
	//=================================================================================================//
	void LocalDirectionalDiffusion::initializeFiberDirection()
	{
		base_particles_->registerVariable(local_bias_direction_, "Fiber");
		base_particles_->addVariableNameToList<Vecd>(reload_local_parameters_, "Fiber");
	}
	//=================================================================================================//
	void LocalDirectionalDiffusion::readFromXmlForLocalParameters(const std::string &filefullpath)
	{
		BaseMaterial::readFromXmlForLocalParameters(filefullpath);
		size_t total_real_particles = base_particles_->total_real_particles_;
		for (size_t i = 0; i != total_real_particles; i++)
		{
			Matd diff_i = diff_cf_ * Matd(1.0) + bias_diff_cf_ * SimTK::outer(local_bias_direction_[i], local_bias_direction_[i]);
			local_transformed_diffusivity_.push_back(inverseCholeskyDecomposition(diff_i));
		}
		std::cout << "\n Local diffusion parameters setup finished " << std::endl;
	};
	//=================================================================================================//
}
