/**
 * @file 	constraint_dynamics.cpp
 * @author	Luhui Han, Chi Zhang and Xiangyu Hu
 */

#include "constrain.h"
#include "vector_functions.h"
#include <numeric>


namespace SPH
{
	//==========================================================================================================//
		Constrain3DSolidBodyRotation::
			Constrain3DSolidBodyRotation(SPHBody &sph_body, Vecd mass_center, Matd inertia_tensor)
			: LocalDynamics(sph_body), SolidDataSimple(sph_body),
			vel_(particles_->vel_), pos_(particles_->pos_), compute_total_moment_of_momentum_(sph_body, mass_center),
			mass_center_(mass_center), moment_of_inertia_(inertia_tensor){}
		//=================================================================================================//
		void Constrain3DSolidBodyRotation::setupDynamics(Real dt)
		{
			angular_velocity_ = moment_of_inertia_.inverse()*compute_total_moment_of_momentum_.parallel_exec(dt);
		}
		//=================================================================================================//
		void Constrain3DSolidBodyRotation::update(size_t index_i, Real dt)
		{
			linear_velocity_ = angular_velocity_.cross( (pos_[index_i] - mass_center_));
			vel_[index_i] -= linear_velocity_;
		}
		//=================================================================================================//
}
