#include "fluid_boundary_static_confinement.h"

namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		StaticConfinementTransportVelocity::StaticConfinementTransportVelocity(NearShapeSurface& near_surface, Real coefficient)
			: BaseLocalDynamics<BodyPartByCell>(near_surface), FluidDataSimple(sph_body_),
			pos_(particles_->pos_), surface_indicator_(particles_->surface_indicator_),
			smoothing_length_sqr_(pow(sph_body_.sph_adaptation_->ReferenceSmoothingLength(), 2)),
			coefficient_(coefficient),
			level_set_shape_(&near_surface.level_set_shape_) {}
		//=================================================================================================//
        void StaticConfinementTransportVelocity::interaction(size_t index_i, Real dt)
		{
			Vecd acceleration_trans = Vecd::Zero();
            /*below for debuging*/
            //Vecd pos_tem = pos_[index_i];

			// acceleration for transport velocity
			acceleration_trans -= 2.0 * level_set_shape_->computeKernelGradientIntegral(pos_[index_i]);
			/** correcting particle position */
            //if (surface_indicator_[index_i] == 1 || surface_indicator_[index_i]==0)
            if (surface_indicator_[index_i] == 0)
            {
                pos_[index_i] += coefficient_ * smoothing_length_sqr_ * acceleration_trans;
            }
			/*std::string output_folder = "./output";
			std::string filefullpath = output_folder + "/" + "transportVelocity_wall_level_" + std::to_string(dt) + ".dat";
			std::ofstream out_file(filefullpath.c_str(), std::ios::app);
			out_file << index_i << " " << surface_indicator_[index_i] << " " << acceleration_trans[0] << " " << acceleration_trans[1]
						<< " " << pos_tem[0] << " " << pos_tem[1] << " " << pos_[index_i][0] << " " << pos_[index_i][1] << " " << dt << std::endl;*/
		}
		//=================================================================================================//
		StaticConfinementViscousAcceleration::StaticConfinementViscousAcceleration(NearShapeSurface& near_surface)
			: BaseLocalDynamics<BodyPartByCell>(near_surface), FluidDataSimple(sph_body_),
			pos_(particles_->pos_), acc_prior_(particles_->acc_prior_), rho_(particles_->rho_),
			mu_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()).ReferenceViscosity()), vel_(particles_->vel_),
			level_set_shape_(&near_surface.level_set_shape_) {}
		//=================================================================================================//
        void StaticConfinementViscousAcceleration::interaction(size_t index_i, Real dt)
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
			StaticConfinementExtendIntegration1stHalf(NearShapeSurface& near_surface, Real penalty_strength)
			: BaseLocalDynamics<BodyPartByCell>(near_surface), FluidDataSimple(sph_body_),
			fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())), 
			rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")),
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
		StaticConfinementIntegration1stHalfPenaltyVelocity::
			StaticConfinementIntegration1stHalfPenaltyVelocity(NearShapeSurface& near_surface, Real sound_speed, Real penalty_strength)
			: BaseLocalDynamics<BodyPartByCell>(near_surface), FluidDataSimple(sph_body_),
			fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())), c_0_(sound_speed),
			rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")),
			pos_(particles_->pos_), vel_(particles_->vel_),
			acc_(particles_->acc_),
			level_set_shape_(&near_surface.level_set_shape_),
			riemann_solver_(fluid_, fluid_), penalty_strength_(penalty_strength) {}
		//=================================================================================================//
		void StaticConfinementIntegration1stHalfPenaltyVelocity::update(size_t index_i, Real dt)
		{
			Vecd kernel_gradient = level_set_shape_->computeKernelGradientIntegral(pos_[index_i]);
			acc_[index_i] -= 2.0 * p_[index_i] * kernel_gradient / rho_[index_i];

			//Real penalty_pressure = 0.5 * c_0_ * c_0_ * rho_[index_i];
			Real penalty_pressure = 0.5 * vel_[index_i].squaredNorm() * rho_[index_i];
			Real distance_to_the_wall = abs(level_set_shape_->findSignedDistance(pos_[index_i]));
			Real ratio = distance_to_the_wall / (3.0 * sph_body_.sph_adaptation_->ReferenceSpacing());
			Real penalty = ratio < 1.0 ? (1.0 - ratio) * (1.0 - ratio) * 0.5 * penalty_pressure : 0.0;

			acc_[index_i] -= 2.0 * penalty_strength_ * penalty * kernel_gradient / rho_[index_i];
		}
		//=================================================================================================//
		StaticConfinementFreeSurfaceIndication::StaticConfinementFreeSurfaceIndication(NearShapeSurface &near_surface)
		: BaseLocalDynamics<BodyPartByCell>(near_surface), FluidDataSimple(sph_body_),
			pos_(particles_->pos_), surface_indicator_(particles_->surface_indicator_),
			level_set_shape_(&near_surface.level_set_shape_), pos_div_(*particles_->getVariableByName<Real>("DensityChangeRate"))
		{}
		//=================================================================================================//
		void StaticConfinementFreeSurfaceIndication::interaction(size_t index_i, Real dt )
		{
			Real pos_div = 0.0;
			Vecd kernel_gradient = level_set_shape_->computeKernelGradientIntegral(pos_[index_i]);
			pos_div -= kernel_gradient.norm() * abs(level_set_shape_->findSignedDistance(pos_[index_i]));
			pos_div_[index_i] += pos_div;
		}
		//=================================================================================================//
		StaticConfinementBounding::StaticConfinementBounding(NearShapeSurface& near_surface)
			: BaseLocalDynamics<BodyPartByCell>(near_surface), FluidDataSimple(sph_body_),
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
        StaticConfinementWithPenalty::StaticConfinementWithPenalty(NearShapeSurface &near_surface, Real sound_speed, Real penalty_strength)
			: density_summation_(near_surface), pressure_relaxation_(near_surface),
			density_relaxation_(near_surface), transport_velocity_(near_surface),
			viscous_acceleration_(near_surface), extend_intergration_1st_half_(near_surface, penalty_strength),
			surface_bounding_(near_surface), extend_intergration_1st_half_Velocity(near_surface, sound_speed, penalty_strength)
		{}
        //=================================================================================================//
        StaticConfinementGeneral::StaticConfinementGeneral(NearShapeSurface &near_surface)
            : density_summation_(near_surface), pressure_relaxation_(near_surface),
                density_relaxation_(near_surface), transport_velocity_(near_surface),
                viscous_acceleration_(near_surface), surface_bounding_(near_surface), free_surface_indication_(near_surface)
        {}
	}
	//=================================================================================================//
}
//=================================================================================================//