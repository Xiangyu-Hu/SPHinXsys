#include "level_set_confinement.h"

namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		StaticConfinementTransportVelocity::StaticConfinementTransportVelocity(NearShapeSurface& near_surface, Real coefficient)
			: BaseLocalDynamics<BodyPartByCell>(near_surface), FluidDataSimple(sph_body_),
			pos_(particles_->pos_), surface_indicator_(particles_->indicator_),
			coefficient_(coefficient), smoothing_length_sqr_(pow(sph_body_.sph_adaptation_->ReferenceSmoothingLength(), 2)),
			level_set_shape_(&near_surface.level_set_shape_), transport_acc_(*particles_->getVariableByName<Vecd>("TransportAcceleration"))
		{}
		//=================================================================================================//
        void StaticConfinementTransportVelocity::update(size_t index_i, Real dt)
		{
			//Vecd acceleration_trans = Vecd::Zero();
            /*below for debuging*/
            //Vecd pos_tem = pos_[index_i];
			Real inv_h_ratio = 1.0;
			// acceleration for transport velocity
			transport_acc_[index_i] -= 2.0 * level_set_shape_->computeKernelGradientIntegral(pos_[index_i]);
			/** correcting particle position */
            //if (surface_indicator_[index_i] == 1 || surface_indicator_[index_i]==0)
            
           /* if (surface_indicator_[index_i] == 0)
            {
                pos_[index_i] += coefficient_ *  smoothing_length_sqr_ * acceleration_trans;
            }*/
            
			 //std::string output_folder = "./output";
				//std::string filefullpath = output_folder + "/" + "transportVelocity_wall_levelset_" + std::to_string(dt) + ".dat";
				//std::ofstream out_file(filefullpath.c_str(), std::ios::app);
				//out_file << pos_[index_i][0] << " " << pos_[index_i][1] << " "<< index_i << " "  << transport_acc_[index_i][0] << " " << transport_acc_[index_i][1]<<" "  << transport_acc_[index_i].norm() << std::endl;
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
			pos_(particles_->pos_), mass_(particles_->mass_), force_prior_(particles_->force_prior_), rho_(particles_->rho_),
			mu_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()).ReferenceViscosity()), vel_(particles_->vel_),
            level_set_shape_(&near_surface.level_set_shape_)
			//, force_from_fluid_(*this->particles_->template registerSharedVariable<Vecd>("ViscousForceFromWall")),
			//kernel_gradient_rij_(*this->particles_->template registerSharedVariable<Real>("KernelGradientRij"))
		{
            particles_->registerVariable(force_from_fluid_, "ViscousForceFromWall"); 
			particles_->registerVariable(kernel_gradient_rij_, "KernelGradientRij");
		}
		//=================================================================================================//
        void StaticConfinementViscousAcceleration::update(size_t index_i, Real dt)
		{
			Vecd force = Vecd::Zero();
			Vecd vel_derivative = Vecd::Zero();
			Vecd vel_level_set_cell_j = Vecd::Zero();
			Real rho_i = rho_[index_i];
			/*Here we give the Level-set boundary velocity as zero, but later we need a vector to set the velocity of each level-set cell*/
			Real phi_r_ij = abs(level_set_shape_->findSignedDistance(pos_[index_i]));
			vel_derivative = 2.0 * (vel_[index_i] - vel_level_set_cell_j);
			Real kernel_gradient_divide_Rij = level_set_shape_->computeKernelGradientDivideRijIntegral(pos_[index_i]);
			force = 2.0 * mu_ * mass_[index_i] * kernel_gradient_divide_Rij * vel_derivative /rho_i;
			force_prior_[index_i] += force;
                        /*below for debuging*/
            force_from_fluid_[index_i] = force;
            kernel_gradient_rij_[index_i] = kernel_gradient_divide_Rij;
			/*for debuging*/
			/*Vecd force = Vecd::Zero();
			force = 2.0 * mu_ * kernel_gradient_divide_Rij * vel_derivative;*/

			std::string output_folder = "./output";
			std::string filefullpath = output_folder + "/" + "viscous_acceleration_wall_levelset_" + std::to_string(dt) + ".dat";
			std::ofstream out_file(filefullpath.c_str(), std::ios::app);
			out_file << this->pos_[index_i][0] << " " << this->pos_[index_i][1] << " "<< index_i << " "  << force[0] << " " << force[1]<<" "  << force.norm() << " "
			<<kernel_gradient_divide_Rij<< std::endl;

		}
		//=================================================================================================//
		BaseForceFromFluidStaticConfinement::BaseForceFromFluidStaticConfinement(NearShapeSurface& near_surface)
			: BaseLocalDynamics<BodyPartByCell>(near_surface), FluidDataSimple(sph_body_), 
			level_set_shape_(&near_surface.level_set_shape_), force_from_fluid_(*this->particles_->template registerSharedVariable<Vecd>("ViscousForceFromFluid"))
		{
		}
		//=================================================================================================//
		ViscousForceFromFluidStaticConfinement::ViscousForceFromFluidStaticConfinement(NearShapeSurface& near_surface)
			:BaseForceFromFluidStaticConfinement(near_surface), pos_(particles_->pos_), rho_(particles_->rho_), mass_(particles_->mass_),
			mu_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()).ReferenceViscosity()), vel_(particles_->vel_)
		{
			//particles_->registerVariable(force_from_fluid_, "ViscousForceFromFluid");
		}
		//=================================================================================================//
		StaticConfinementExtendIntegration1stHalf::
			StaticConfinementExtendIntegration1stHalf(NearShapeSurface& near_surface, Real penalty_strength)
			: BaseLocalDynamics<BodyPartByCell>(near_surface), FluidDataSimple(sph_body_),
			fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())), 
			rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")),
			pos_(particles_->pos_), vel_(particles_->vel_), mass_(particles_->mass_),
			force_(particles_->force_),
			level_set_shape_(&near_surface.level_set_shape_),
			riemann_solver_(fluid_, fluid_), penalty_strength_(penalty_strength) {}
		//=================================================================================================//
		void StaticConfinementExtendIntegration1stHalf::update(size_t index_i, Real dt)
		{
			Vecd kernel_gradient = level_set_shape_->computeKernelGradientIntegral(pos_[index_i]);
			force_[index_i] -= 2.0 * p_[index_i] * kernel_gradient / rho_[index_i];
			
			Real penalty_pressure = p_[index_i];
			Real distance_to_the_wall = abs(level_set_shape_->findSignedDistance(pos_[index_i]));
			Real ratio = distance_to_the_wall  / (3.0 * sph_body_.sph_adaptation_->ReferenceSpacing());
			Real penalty = ratio < 1.0 ? (1.0 - ratio) * (1.0 - ratio) * 1.0 * penalty_pressure: 0.0;
			
			force_[index_i] -= 2.0 * penalty_strength_* penalty * kernel_gradient / rho_[index_i];
		}
		//=================================================================================================//
		StaticConfinementIntegration1stHalfPenaltyVelocity::
			StaticConfinementIntegration1stHalfPenaltyVelocity(NearShapeSurface& near_surface, Real sound_speed, Real penalty_strength)
			: BaseLocalDynamics<BodyPartByCell>(near_surface), FluidDataSimple(sph_body_),
			fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())), c_0_(sound_speed),
			rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")),
			pos_(particles_->pos_), vel_(particles_->vel_), mass_(particles_->mass_),
			force_(particles_->force_),
			level_set_shape_(&near_surface.level_set_shape_),
			riemann_solver_(fluid_, fluid_), penalty_strength_(penalty_strength) {}
		//=================================================================================================//
		void StaticConfinementIntegration1stHalfPenaltyVelocity::update(size_t index_i, Real dt)
		{
			Vecd kernel_gradient = level_set_shape_->computeKernelGradientIntegral(pos_[index_i]);
			force_[index_i] -= 2.0 * p_[index_i] * kernel_gradient / rho_[index_i];

			//Real penalty_pressure = 0.5 * c_0_ * c_0_ * rho_[index_i];
			Real penalty_pressure = 0.5 * vel_[index_i].squaredNorm() * rho_[index_i];
			Real distance_to_the_wall = abs(level_set_shape_->findSignedDistance(pos_[index_i]));
			Real ratio = distance_to_the_wall / (3.0 * sph_body_.sph_adaptation_->ReferenceSpacing());
			Real penalty = ratio < 1.0 ? (1.0 - ratio) * (1.0 - ratio) * 0.5 * penalty_pressure : 0.0;

			force_[index_i] -= 2.0 * penalty_strength_ * penalty * kernel_gradient / rho_[index_i];
		}
		//=================================================================================================//
		StaticConfinementFreeSurfaceIndication::StaticConfinementFreeSurfaceIndication(NearShapeSurface &near_surface)
		: BaseLocalDynamics<BodyPartByCell>(near_surface), FluidDataSimple(sph_body_),
			pos_(particles_->pos_), surface_indicator_(particles_->indicator_),
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
		
	}
	//=================================================================================================//
	
	//=================================================================================================//
}
//=================================================================================================//