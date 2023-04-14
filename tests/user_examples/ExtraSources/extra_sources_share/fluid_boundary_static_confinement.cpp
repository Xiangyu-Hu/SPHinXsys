#include "fluid_boundary_static_confinement.h"

namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		StaticConfinementTransportVelocity::StaticConfinementTransportVelocity(NearShapeSurface& near_surface, Real coefficient)
			: LocalDynamics(near_surface.getSPHBody()), FluidDataSimple(sph_body_),
			pos_(particles_->pos_), surface_indicator_(particles_->surface_indicator_),
			smoothing_length_sqr_(pow(sph_body_.sph_adaptation_->ReferenceSmoothingLength(), 2)),
			coefficient_(coefficient),
			level_set_shape_(&near_surface.level_set_shape_) {}
		//=================================================================================================//
		void StaticConfinementTransportVelocity::update(size_t index_i, Real dt)
		{
			Vecd acceleration_trans = Vecd::Zero();
			// acceleration for transport velocity
			acceleration_trans -= 2.0 * level_set_shape_->computeKernelGradientIntegral(pos_[index_i]);
			/** correcting particle position */
			if (surface_indicator_[index_i] == 0)
				pos_[index_i] += coefficient_ * smoothing_length_sqr_ * acceleration_trans;
		}
		//=================================================================================================//
		StaticConfinementViscousAcceleration::StaticConfinementViscousAcceleration(NearShapeSurface& near_surface)
			: LocalDynamics(near_surface.getSPHBody()), FluidDataSimple(sph_body_),
			pos_(particles_->pos_), acc_prior_(particles_->acc_prior_), rho_(particles_->rho_),
			mu_(particles_->fluid_.ReferenceViscosity()), vel_(particles_->vel_),
			level_set_shape_(&near_surface.level_set_shape_) {}
		//=================================================================================================//
		void StaticConfinementViscousAcceleration::update(size_t index_i, Real dt)
		{
			Vecd acceleration = Vecd::Zero();
			Vecd vel_derivative = Vecd::Zero();
			Vecd vel_level_set_cell_j = Vecd::Zero();
			/*Here we give the Level-set boundary velocity as zero, but later we need a vector to set the velocity of each level-set cell*/
			Real phi_r_ij = abs(level_set_shape_->findSignedDistance(pos_[index_i]));
			vel_derivative = (vel_[index_i] - vel_level_set_cell_j) / (phi_r_ij + TinyReal);
			Vecd kernel_gradient = level_set_shape_->computeKernelGradientIntegral(pos_[index_i]);
			acceleration += 2.0 * mu_ * kernel_gradient.norm() * vel_derivative;
			acc_prior_[index_i] += acceleration / rho_[index_i];
		}
		//=================================================================================================//
		StaticConfinementExtendIntegration1stHalf::
			StaticConfinementExtendIntegration1stHalf(NearShapeSurface& near_surface, Real  sound_speed, Real penalty_strength)
			: LocalDynamics(near_surface.getSPHBody()), FluidDataSimple(sph_body_),
			fluid_(particles_->fluid_), c_0_ (sound_speed),
			rho_(particles_->rho_), p_(particles_->p_),
			pos_(particles_->pos_), vel_(particles_->vel_),
			acc_(particles_->acc_),
			level_set_shape_(&near_surface.level_set_shape_),
			riemann_solver_(fluid_, fluid_), penalty_strength_(penalty_strength) {}
		//=================================================================================================//
		void StaticConfinementExtendIntegration1stHalf::update(size_t index_i, Real dt)
		{
			Vecd kernel_gradient = level_set_shape_->computeKernelGradientIntegral(pos_[index_i]);
			acc_[index_i] -= 2.0 * p_[index_i] * kernel_gradient / rho_[index_i];
			
			Real penalty_pressure = p_[index_i];
			Real distance_to_the_wall = abs(level_set_shape_->findSignedDistance(pos_[index_i]));
			Real ratio = distance_to_the_wall  / (3.0 * sph_body_.sph_adaptation_->ReferenceSpacing());
			Real penalty = ratio < 1.0 ? (1.0 - ratio) * (1.0 - ratio) * 1.0 * penalty_pressure: 0.0;
			
			acc_[index_i] -= 2.0 * penalty_strength_* penalty * kernel_gradient / rho_[index_i];
		}
		//=================================================================================================//
		StaticConfinementBounding::StaticConfinementBounding(NearShapeSurface& near_surface)
			: LocalDynamics(near_surface.getSPHBody()), FluidDataSimple(sph_body_),
			pos_(particles_->pos_),
			constrained_distance_(0.5 * sph_body_.sph_adaptation_->MinimumSpacing())
		{
			level_set_shape_ = &near_surface.level_set_shape_;
		}
		//=================================================================================================//
		void StaticConfinementBounding::update(size_t index_i, Real dt)
		{
			Real phi = level_set_shape_->findSignedDistance(pos_[index_i]);

			if (phi > -constrained_distance_)
			{
				Vecd unit_normal = level_set_shape_->findNormalDirection(pos_[index_i]);
				pos_[index_i] -= (phi + constrained_distance_) * unit_normal;
			}
		}
		//=================================================================================================//
		StaticConfinementWithBounding::StaticConfinementWithBounding(NearShapeSurface& near_surface)
			: density_summation_(near_surface), pressure_relaxation_(near_surface),
			density_relaxation_(near_surface), surface_bounding_(near_surface)
		{}
		//=================================================================================================//
		StaticConfinementWithPenalty::StaticConfinementWithPenalty(NearShapeSurface& near_surface, Real sound_speed, Real penalty_strength)
			: density_summation_(near_surface), pressure_relaxation_(near_surface),
			density_relaxation_(near_surface), transport_velocity_(near_surface),
			viscous_acceleration_(near_surface), extend_intergration_1st_half_(near_surface, sound_speed, penalty_strength),
			surface_bounding_(near_surface)
		{}
	}
	//=================================================================================================//
}
//=================================================================================================//