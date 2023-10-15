#include "level_set_confinement.h"

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
			 //std::string output_folder = "./output";
				//std::string filefullpath = output_folder + "/" + "transportVelocity_wall_levelset_" + std::to_string(dt) + ".dat";
				//std::ofstream out_file(filefullpath.c_str(), std::ios::app);
				//out_file << pos_[index_i][0] << " " << pos_[index_i][1] << " "<< index_i << " "  << acceleration_trans[0] << " " << acceleration_trans[1]<<" "  << acceleration_trans.norm() << std::endl;
				///** correcting particle position */

			/*std::string output_folder = "./output";
			std::string filefullpath = output_folder + "/" + "transportVelocity_wall_levelset_" + std::to_string(dt) + ".dat";
			std::ofstream out_file(filefullpath.c_str(), std::ios::app);
			out_file <<this->particles_->pos_[index_i][0]<<" "<<this->particles_->pos_[index_i][1]<<" " <<index_i<< "  "<<  acceleration_trans.norm()<<std::endl;
			out_file << " \n";*/
		}
		//=================================================================================================//
		StaticConfinementViscousAcceleration::StaticConfinementViscousAcceleration(NearShapeSurface& near_surface)
			: BaseLocalDynamics<BodyPartByCell>(near_surface), FluidDataSimple(sph_body_),
			pos_(particles_->pos_), acc_prior_(particles_->acc_prior_), rho_(particles_->rho_),
			mu_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()).ReferenceViscosity()), vel_(particles_->vel_),
			level_set_shape_(&near_surface.level_set_shape_) {}
		//=================================================================================================//
        void StaticConfinementViscousAcceleration::update(size_t index_i, Real dt)
		{
			Vecd acceleration = Vecd::Zero();
			Vecd vel_derivative = Vecd::Zero();
			Vecd vel_level_set_cell_j = Vecd::Zero();
			Real rho_i = rho_[index_i];
			/*Here we give the Level-set boundary velocity as zero, but later we need a vector to set the velocity of each level-set cell*/
			Real phi_r_ij = abs(level_set_shape_->findSignedDistance(pos_[index_i]));
			vel_derivative = 2.0 * (vel_[index_i] - vel_level_set_cell_j);
			Real kernel_gradient_divide_Rij = level_set_shape_->computeKernelGradientDivideRijIntegral(pos_[index_i]);
			acceleration += 2.0 * mu_ * kernel_gradient_divide_Rij * vel_derivative /rho_i;
			acc_prior_[index_i] += acceleration / rho_[index_i];

				/*std::string output_folder = "./output";
				std::string filefullpath = output_folder + "/" + "viscous_acceleration_wall_levelset_" + std::to_string(dt) + ".dat";
				std::ofstream out_file(filefullpath.c_str(), std::ios::app);
				out_file << this->pos_[index_i][0] << " " << this->pos_[index_i][1] << " "<< index_i << " "  << acceleration[0] << " " << acceleration[1]<<" "  << acceleration.norm() << " "<<kernel_gradient_divide_Rij<< std::endl;*/
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
			level_set_shape_(&near_surface.level_set_shape_), pos_div_(*particles_->getVariableByName<Real>("PositionDivergence"))
		{}
		//=================================================================================================//
		void StaticConfinementFreeSurfaceIndication::interaction(size_t index_i, Real dt )
		{
			//std::string output_folder_1 = "./output";
			//std::string filefullpath_1 = output_folder_1 + "/" + "position_divergence_levelset_before" + std::to_string(dt) + ".dat";
			//std::ofstream out_file_1(filefullpath_1.c_str(), std::ios::app);
			//out_file_1 <<pos_[index_i][0]<<" "<<pos_[index_i][1]<<" " <<index_i<< "  "<<  pos_div_[index_i]<<std::endl;
			////out_file << std::fixed << std::setprecision(2)<<pos_[index_i][0]<< "  " <<pos_[index_i][1]<< "  "<< pos_div << "  "<<  pos_div_[index_i]<<std::endl;
			//out_file_1 << " \n";

			Real pos_div = - level_set_shape_->computeKernelGradientMultiplyRijIntegral(pos_[index_i]);
			pos_div_[index_i] += pos_div;
			//std::string output_folder = "./output";
			//std::string filefullpath = output_folder + "/" + "position_divergence_levelset_after" + std::to_string(dt) + ".dat";
			//std::ofstream out_file(filefullpath.c_str(), std::ios::app);
			//out_file <<pos_[index_i][0]<<" "<<pos_[index_i][1]<<" " <<index_i<< "  "<<  pos_div_[index_i]<<std::endl;
			////out_file << std::fixed << std::setprecision(2)<<pos_[index_i][0]<< "  " <<pos_[index_i][1]<< "  "<< pos_div << "  "<<  pos_div_[index_i]<<std::endl;
			//out_file << " \n";
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
		//=================================================================================================//
		MovingConfinementDensitySummation::MovingConfinementDensitySummation(NearShapeSurfaceTracing& near_surface_tracing )
			: BaseLocalDynamics<BodyPartByCell>(near_surface_tracing), FluidDataSimple(sph_body_),
		rho0_(sph_body_.base_material_->ReferenceDensity()), near_surface_tracing_(near_surface_tracing),
		inv_sigma0_(1.0 / sph_body_.sph_adaptation_->LatticeNumberDensity()),
		mass_(particles_->mass_), rho_sum_(*particles_->getVariableByName<Real>("DensitySummation")), pos_(particles_->pos_),
		level_set_shape_(&near_surface_tracing.level_set_shape_){}
		//=================================================================================================//
		void MovingConfinementDensitySummation::update(size_t index_i, Real dt)
		{
			Real inv_Vol_0_i = rho0_ / mass_[index_i];
			rho_sum_[index_i] +=
				level_set_shape_->computeKernelIntegral(near_surface_tracing_.tracing_cell_method_base_.tracingPosition(pos_[index_i])) * inv_Vol_0_i * rho0_ * inv_sigma0_;
		}
		//=================================================================================================//
		MovingConfinementIntegration1stHalf::MovingConfinementIntegration1stHalf(NearShapeSurfaceTracing& near_surface_tracing)
			 : BaseLocalDynamics<BodyPartByCell>(near_surface_tracing), FluidDataSimple(sph_body_),
		fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())), near_surface_tracing_(near_surface_tracing),
		rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")),
		pos_(particles_->pos_), vel_(particles_->vel_),
		acc_(particles_->acc_),
		level_set_shape_(&near_surface_tracing.level_set_shape_),
		riemann_solver_(fluid_, fluid_) {}
		//=================================================================================================//
		void MovingConfinementIntegration1stHalf::update(size_t index_i, Real dt)
		{
			Vecd kernel_gradient_previous = level_set_shape_->computeKernelGradientIntegral(near_surface_tracing_.tracing_cell_method_base_.tracingPosition(pos_[index_i]));
			Vecd kernel_gradient = near_surface_tracing_.tracing_cell_method_base_.updateNormalForVector(kernel_gradient_previous);
			acc_[index_i] -= 2.0 * p_[index_i] * kernel_gradient / rho_[index_i];
		}
		//=================================================================================================//
		MovingConfinementIntegration2ndHalf::MovingConfinementIntegration2ndHalf(NearShapeSurfaceTracing& near_surface_tracing)
			 : BaseLocalDynamics<BodyPartByCell>(near_surface_tracing), FluidDataSimple(sph_body_),
		fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
		rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")), drho_dt_(*particles_->getVariableByName<Real>("DensityChangeRate")),
		pos_(particles_->pos_), vel_(particles_->vel_), near_surface_tracing_(near_surface_tracing),
		level_set_shape_(&near_surface_tracing.level_set_shape_),
		riemann_solver_(fluid_, fluid_) {}
		//=================================================================================================//
		void MovingConfinementIntegration2ndHalf::update(size_t index_i, Real dt)
		{
			Vecd kernel_gradient_previous = level_set_shape_->computeKernelGradientIntegral(near_surface_tracing_.tracing_cell_method_base_.tracingPosition(pos_[index_i]));
			Vecd kernel_gradient = near_surface_tracing_.tracing_cell_method_base_.updateNormalForVector(kernel_gradient_previous);
			Vecd vel_in_wall = -vel_[index_i];
			drho_dt_[index_i] += rho_[index_i] * (vel_[index_i] - vel_in_wall).dot(kernel_gradient);
		}
		//=================================================================================================//
		//template <class TracingMethodType>
		MovingConfinementBounding::MovingConfinementBounding(NearShapeSurfaceTracing& near_surface_tracing)
			:BaseLocalDynamics<BodyPartByCell>(near_surface_tracing), FluidDataSimple(sph_body_),
			pos_(particles_->pos_),level_set_shape_(&near_surface_tracing.level_set_shape_), 
			near_surface_tracing_(near_surface_tracing),
			constrained_distance_(0.5 * sph_body_.sph_adaptation_->MinimumSpacing()){}
		//=================================================================================================//
		void MovingConfinementBounding::update(size_t index_i, Real dt)
		{
			Real phi = level_set_shape_->findSignedDistance(near_surface_tracing_.tracing_cell_method_base_.tracingPosition(pos_[index_i]));

			if (phi > -constrained_distance_)
			{
				Vecd unit_normal_previous = level_set_shape_->findNormalDirection(near_surface_tracing_.tracing_cell_method_base_.tracingPosition(pos_[index_i]));
				Vecd unit_normal = near_surface_tracing_.tracing_cell_method_base_.updateNormalForVector(unit_normal_previous);
				pos_[index_i] -= (phi + constrained_distance_) * unit_normal;
			}
		}
		//=================================================================================================//
		MovingConfinementFreeSurfaceIndication::MovingConfinementFreeSurfaceIndication(NearShapeSurfaceTracing& near_surface_tracing)
		: BaseLocalDynamics<BodyPartByCell>(near_surface_tracing), FluidDataSimple(sph_body_),
			pos_(particles_->pos_), surface_indicator_(particles_->surface_indicator_),
			near_surface_tracing_(near_surface_tracing),
			level_set_shape_(&near_surface_tracing.level_set_shape_), pos_div_(*particles_->getVariableByName<Real>("DensityChangeRate"))
		{}
		//=================================================================================================//
		void MovingConfinementFreeSurfaceIndication::interaction(size_t index_i, Real dt )
		{
			Real pos_div = - level_set_shape_->computeKernelGradientMultiplyRijIntegral(near_surface_tracing_.tracing_cell_method_base_.tracingPosition(pos_[index_i]));
			pos_div_[index_i] += pos_div;
		}
		//=================================================================================================//
		 MovingConfinementTransportVelocity::MovingConfinementTransportVelocity(NearShapeSurfaceTracing& near_surface_tracing, Real coefficient)
			: BaseLocalDynamics<BodyPartByCell>(near_surface_tracing), FluidDataSimple(sph_body_),
			pos_(particles_->pos_), surface_indicator_(particles_->surface_indicator_),
			smoothing_length_sqr_(pow(sph_body_.sph_adaptation_->ReferenceSmoothingLength(), 2)),
			coefficient_(coefficient), near_surface_tracing_(near_surface_tracing),
			level_set_shape_(&near_surface_tracing.level_set_shape_) {}
		//=================================================================================================//
        void MovingConfinementTransportVelocity::interaction(size_t index_i, Real dt)
		{
			Vecd acceleration_trans = Vecd::Zero();
			Vecd kernel_gradient_previous = level_set_shape_->computeKernelGradientIntegral(near_surface_tracing_.tracing_cell_method_base_.tracingPosition(pos_[index_i]));
			Vecd kernel_gradient = near_surface_tracing_.tracing_cell_method_base_.updateNormalForVector(kernel_gradient_previous);
			acceleration_trans -= 2.0 * kernel_gradient;
            if (surface_indicator_[index_i] == 0)
            {
                pos_[index_i] += coefficient_ * smoothing_length_sqr_ * acceleration_trans;
            }
		}
		//=================================================================================================//
		MovingConfinementViscousAcceleration::MovingConfinementViscousAcceleration(NearShapeSurfaceTracing& near_surface_tracing)
			:BaseLocalDynamics<BodyPartByCell>(near_surface_tracing), FluidDataSimple(sph_body_),
			pos_(particles_->pos_), acc_prior_(particles_->acc_prior_), rho_(particles_->rho_),
			mu_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()).ReferenceViscosity()), vel_(particles_->vel_),
			near_surface_tracing_(near_surface_tracing),
			level_set_shape_(&near_surface_tracing.level_set_shape_) {}
		//=================================================================================================//
		void MovingConfinementViscousAcceleration::update(size_t index_i, Real dt)
		{
			Vecd acceleration = Vecd::Zero();
			Vecd vel_derivative = Vecd::Zero();
			Vecd vel_level_set_cell_j = Vecd::Zero();
			Real rho_i = rho_[index_i];
			/*Here we give the Level-set boundary velocity as zero, but later we need a vector to set the velocity of each level-set cell*/
			Real phi_r_ij = abs(level_set_shape_->findSignedDistance(pos_[index_i]));
			vel_derivative = 2.0 * (vel_[index_i] - vel_level_set_cell_j);
			Real kernel_gradient_divide_Rij = level_set_shape_->computeKernelGradientDivideRijIntegral(near_surface_tracing_.tracing_cell_method_base_.tracingPosition(pos_[index_i]));
			acceleration += 2.0 * mu_ * kernel_gradient_divide_Rij * vel_derivative /rho_i;
			acc_prior_[index_i] += acceleration / rho_[index_i];

		}
		//=================================================================================================//
		MovingConfinementGeneral::MovingConfinementGeneral(NearShapeSurfaceTracing& near_surface_tracing)
		: density_summation_(near_surface_tracing), pressure_relaxation_(near_surface_tracing),
		  density_relaxation_(near_surface_tracing), surface_bounding_(near_surface_tracing), 
	      transport_velocity_(near_surface_tracing),free_surface_indication_(near_surface_tracing),
			viscous_acceleration_(near_surface_tracing)
		{}
	}
	//=================================================================================================//
	
	//=================================================================================================//
}
//=================================================================================================//