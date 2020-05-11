/**
 * @file 	active_muscle_dynamics.cpp
 * @brief 	In is file, we define functions decleared in active muscle dynamics.h
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.3.1
 *			Here, we need identify the physical differences between electrophysilogy and active muscle.
 *			The former is on the diffusion and electro-chemical reaction happens in tissue.
 *			The latter is for muscle dynamics which is driven by an external injection of energy.
 *			As the naming of class, function and variables in this code need be based on physics.
 *			We will identify physical differences by properly choosing names.
 *			Xiangyu Hu
 */
#include "active_muscle_dynamics.h"

using namespace SimTK;

namespace SPH
{
	namespace active_muscle_dynamics
	{
		//=================================================================================================//
		Vecd SpringConstrainMuscleRegion::GetAcceleration(Vecd &disp, Real mass)
		{
			Vecd spring_force(0);
			for(int i = 0; i < disp.size(); i++)
			{
				spring_force[i] = -stiffness_[i] * disp[i] / mass;
			}
			return spring_force;
		}
		//=================================================================================================//
		void SpringConstrainMuscleRegion::ConstraintAParticle(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i
				= particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i 
				= particles_->elastic_body_data_[index_particle_i];

			Vecd disp_from_0 = base_particle_data_i.pos_n_ - base_particle_data_i.pos_0_;
			base_particle_data_i.vel_n_ 	+=  dt * GetAcceleration(disp_from_0, elastic_data_i.mass_);
			base_particle_data_i.pos_n_ 	+=  dt * dt * GetAcceleration(disp_from_0, elastic_data_i.mass_);
		}
		//=================================================================================================//
		void ImposingStress
			::ConstraintAParticle(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];
			ActiveMuscleData &active_muscle_data_i = particles_->active_muscle_data_[index_particle_i];

			active_muscle_data_i.active_stress_= getStress(base_particle_data_i.pos_0_);
		}
		//=================================================================================================//
    }
}
