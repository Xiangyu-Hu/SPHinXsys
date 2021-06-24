/**
* @file 	complex_solid.hpp
* @brief 	These are implementation of the classes for complex solids.
* @author	Xiangyu Hu and Chi Zhang
*/
#pragma once

#include "complex_solid.h"
#include "solid_particles.h"

using namespace std;

namespace SPH {
	//=============================================================================================//
	template<class MuscleType>
    ActiveMuscle<MuscleType>::ActiveMuscle() : MuscleType(), active_muscle_particles_(NULL) 
	{
		MuscleType::material_name_ = "ActiveMuscle";
	}
 	//=============================================================================================//
	template<class MuscleType>
    void ActiveMuscle<MuscleType>::assignDerivedMaterialParameters() 
	{
		MuscleType::assignDerivedMaterialParameters();
	}
  	//=============================================================================================//
	template<class MuscleType>
    void ActiveMuscle<MuscleType>::
        assignActiveMuscleParticles(ActiveMuscleParticles* active_muscle_particles)
	{
		active_muscle_particles_ = active_muscle_particles;
	}
  	//=============================================================================================//
	template<class MuscleType>
    Matd ActiveMuscle<MuscleType>::ConstitutiveRelation(Matd& deformation, size_t index_i)
	{
		return MuscleType::ConstitutiveRelation(deformation, index_i) + 
			active_muscle_particles_->active_contraction_stress_[index_i] * MuscleType::MuscleFiberDirection(index_i);
	}
  //=============================================================================================//
 }
