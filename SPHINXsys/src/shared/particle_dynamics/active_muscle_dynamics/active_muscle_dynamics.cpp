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
		MuscleActivation::
			MuscleActivation(SolidBody* body) :
			ParticleDynamicsSimple(body), ActiveMuscleDataDelegateSimple(body),
			pos_0_(particles_->pos_0_), active_contraction_stress_(particles_->active_contraction_stress_)
		{};
		//=================================================================================================//
		SpringConstrainMuscleRegion::
			SpringConstrainMuscleRegion(SolidBody* body, BodyPartByParticle* body_part) :
			PartDynamicsByParticle(body, body_part),
			ActiveMuscleDataDelegateSimple(body), mass_(particles_->mass_),
			pos_n_(particles_->pos_n_), pos_0_(particles_->pos_0_),
			vel_n_(particles_->vel_n_)
		{
		}
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
		void SpringConstrainMuscleRegion::Update(size_t index_i, Real dt)
		{
			Vecd disp_from_0 = pos_n_[index_i] - pos_0_[index_i];
			vel_n_[index_i] +=  dt * GetAcceleration(disp_from_0, mass_[index_i]);
			pos_n_[index_i] +=  dt * dt * GetAcceleration(disp_from_0, mass_[index_i]);
		}
		//=================================================================================================//
		ImposingStress::
			ImposingStress(SolidBody* body, SolidBodyPartForSimbody* body_part) :
			PartDynamicsByParticle(body, body_part),
			ActiveMuscleDataDelegateSimple(body),
			pos_0_(particles_->pos_0_), active_stress_(particles_->active_stress_)
		{
		}
		//=================================================================================================//
		void ImposingStress
			::Update(size_t index_i, Real dt)
		{
			active_stress_[index_i] = getStress(pos_0_[index_i]);
		}
		//=================================================================================================//
    }
}
