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
			ParticleDynamicsSimple(solid_body), ActiveMuscleDataDelegateSimple(solid_body),
			pos_0_(particles_->pos_0_), active_contraction_stress_(particles_->active_contraction_stress_) {};
		//=================================================================================================//
		SpringConstrainMuscleRegion::
			SpringConstrainMuscleRegion(SolidBody &solid_body, BodyPartByParticle &body_part) :
			PartSimpleDynamicsByParticle(solid_body, body_part),
			ActiveMuscleDataDelegateSimple(solid_body), mass_(particles_->mass_),
			pos_n_(particles_->pos_n_), pos_0_(particles_->pos_0_),
			vel_n_(particles_->vel_n_) {}
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
			Vecd disp_from_0 = pos_n_[index_i] - pos_0_[index_i];
			vel_n_[index_i] +=  dt * getAcceleration(disp_from_0, mass_[index_i]);
			pos_n_[index_i] +=  dt * dt * getAcceleration(disp_from_0, mass_[index_i]);
		}
		//=================================================================================================//
		ImposingStress::
			ImposingStress(SolidBody &solid_body, SolidBodyPartForSimbody &body_part) :
			PartSimpleDynamicsByParticle(solid_body, body_part),
			ActiveMuscleDataDelegateSimple(solid_body),
			pos_0_(particles_->pos_0_), active_stress_(particles_->active_stress_) {}
		//=================================================================================================//
		void ImposingStress
			::Update(size_t index_i, Real dt)
		{
			active_stress_[index_i] = getStress(pos_0_[index_i]);
		}
		//=================================================================================================//
    }
}
