#pragma once

#include "base_particles.hpp"
#include "complex_solid.h"

namespace SPH
{
//=============================================================================================//
template <class MuscleType>
template <typename... ConstructorArgs>
ActiveMuscle<MuscleType>::ActiveMuscle(ConstructorArgs &&...args)
    : MuscleType(std::forward<ConstructorArgs>(args)...)
{
    MuscleType::material_type_name_ = "ActiveMuscle";
}
//=============================================================================================//
template <class MuscleType>
void ActiveMuscle<MuscleType>::initializeLocalParameters(BaseParticles *base_particles)
{
    MuscleType::initializeLocalParameters(base_particles);
    base_particles->registerVariable(active_contraction_stress_, "ActiveContractionStress");
}
//=============================================================================================//
template <class MuscleType>
Matd ActiveMuscle<MuscleType>::StressPK2(Matd &deformation, size_t index_i)
{
    return MuscleType::StressPK2(deformation, index_i) +
           active_contraction_stress_[index_i] * MuscleType::MuscleFiberDirection(index_i);
}
//=============================================================================================//
} // namespace SPH
