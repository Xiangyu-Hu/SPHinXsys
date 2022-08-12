/**
 * @file 	active_muscle_dynamics.cpp
 * @brief 	In is file, we define functions decleared in active muscle dynamics.h
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "active_muscle_dynamics.h"

using namespace SimTK;

namespace SPH
{
	namespace active_muscle_dynamics
	{
		//=================================================================================================//
		MuscleActivation::
			MuscleActivation(SolidBody &solid_body) :
			ParticleDynamicsSimple(solid_body), ElasticSolidDataSimple(solid_body),
			pos0_(particles_->pos0_), 
			active_contraction_stress_(*particles_->getVariableByName<Real>("ActiveContractionStress")) {};
		//=================================================================================================//
		SpringConstrainMuscleRegion::
			SpringConstrainMuscleRegion(SolidBody &solid_body, BodyPartByParticle &body_part) :
			PartSimpleDynamicsByParticle(solid_body, body_part),
			ElasticSolidDataSimple(solid_body), mass_(particles_->mass_),
			pos_(particles_->pos_), pos0_(particles_->pos0_),
			vel_(particles_->vel_) {}
		//=================================================================================================//
		Vecd SpringConstrainMuscleRegion::getAcceleration(Vecd &disp, Real mass)
		{
			Vecd spring_force(0);
			for(int i = 0; i < disp.size(); i++)
			{
				spring_force[i] = -stiffness_[i] * disp[i] / mass;
			}
			return spring_force;
		}
		//=================================================================================================//
		void SpringConstrainMuscleRegion::Update(size_t index_i, Real dt)
		{
			Vecd disp_from_0 = pos_[index_i] - pos0_[index_i];
			vel_[index_i] +=  dt * getAcceleration(disp_from_0, mass_[index_i]);
			pos_[index_i] +=  dt * dt * getAcceleration(disp_from_0, mass_[index_i]);
		}
		//=================================================================================================//
    }
}
