/**
 * @file 	constraint_dynamics.cpp
 * @author	Luhui Han, Chi Zhang and Xiangyu Hu
 */

#include "constrain.h"

#include <numeric>


namespace SPH
{
		Constrain2DSolidBodyRotation::
			Constrain2DSolidBodyRotation(SPHBody &sph_body, Vecd mass_center, Real moment_of_inertia)
			: LocalDynamics(sph_body), SolidDataSimple(sph_body),
			vel_(particles_->vel_), pos_(particles_->pos_), compute_total_moment_of_momentum_(sph_body, mass_center),
			mass_center_(mass_center), moment_of_inertia_(moment_of_inertia) {}
		//=================================================================================================//
		void Constrain2DSolidBodyRotation::setupDynamics(Real dt)
		{
			angular_velocity_ = compute_total_moment_of_momentum_.parallel_exec(dt) / moment_of_inertia_;
		}
		//=================================================================================================//
		void Constrain2DSolidBodyRotation::update(size_t index_i, Real dt)
		{
			Real x = pos_[index_i][0] - mass_center_[0];
			Real y = pos_[index_i][1] - mass_center_[1];
			Real local_radius = sqrt(pow(y, 2.0) + pow(x, 2.0));
			Real angular = atan2(y, x);

			linear_velocity_[1] = angular_velocity_ * local_radius * cos(angular);
			linear_velocity_[0] = -angular_velocity_ * local_radius * sin(angular);
			vel_[index_i] -= linear_velocity_;
		}
		//=================================================================================================//
}
