#pragma once

#include "base_particles.hpp"
#include "complex_solid.h"

namespace SPH
{
//=============================================================================================//
template <class MuscleType>
template <typename... Args>
ActiveMuscle<MuscleType>::ActiveMuscle(Args &&...args)
    : MuscleType(std::forward<Args>(args)...),
      active_contraction_stress_(nullptr)
{
    MuscleType::material_type_name_ = "ActiveMuscle";
}
//=============================================================================================//
template <class MuscleType>
void ActiveMuscle<MuscleType>::initializeLocalParameters(BaseParticles *base_particles)
{
    MuscleType::initializeLocalParameters(base_particles);
    active_contraction_stress_ = base_particles->registerSharedVariable<Real>("ActiveContractionStress");
}
//=============================================================================================//
template <class MuscleType>
Matd ActiveMuscle<MuscleType>::StressPK2(Matd &deformation, size_t index_i)
{
    return MuscleType::StressPK2(deformation, index_i) +
           (*active_contraction_stress_)[index_i] * MuscleType::MuscleFiberDirection(index_i);
}
//=============================================================================================//
} // namespace SPH
