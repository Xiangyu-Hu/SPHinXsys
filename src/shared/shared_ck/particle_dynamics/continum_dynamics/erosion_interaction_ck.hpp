#ifndef EROSION_INTERACTION_CK_HPP
#define EROSION_INTERACTION_CK_HPP

#include "erosion_interaction_ck.h"

namespace SPH
{
//=================================================================================================//
template <class SoilContactRelationType>
Interaction<Soil>::Interaction(SoilContactRelationType &soil_contact_relation)
{
    StdVec<BaseParticles *> contact_particles = soil_contact_relation.getContactParticles();
    for (size_t k = 0; k != contact_particles.size(); ++k)
    {
        dv_soil_n_.push_back(contact_particles[k]->template getVariableByName<Vecd>("NormalDirection"));
        dv_soil_vel_.push_back(contact_particles[k]->template getVariableByName<Vecd>("Velocity"));
        dv_soil_Vol_.push_back(contact_particles[k]->template getVariableByName<Real>("VolumetricMeasure"));
        dv_soil_p_.push_back(contact_particles[k]->template getVariableByName<Real>("Pressure"));
        dv_soil_shear_stress_tensor_.push_back(contact_particles[k]->template getVariableByName<Mat3d>("ShearStressTensor"));
    }
}
//=================================================================================================//
template <class FluidContactRelationType>
Interaction<Fluid>::Interaction(FluidContactRelationType &fluid_contact_relation)
{
    StdVec<BaseParticles *> contact_particles = fluid_contact_relation.getContactParticles();
    for (size_t k = 0; k != contact_particles.size(); ++k)
    {
        dv_fluid_vel_.push_back(contact_particles[k]->template registerStateVariableOnly<Vecd>("Velocity"));
        dv_fluid_Vol_.push_back(contact_particles[k]->template registerStateVariableOnly<Real>("VolumetricMeasure"));
        dv_fluid_p_.push_back(contact_particles[k]->template registerStateVariableOnly<Real>("Pressure"));
    }
}
//=================================================================================================//
} // namespace SPH
#endif // EROSION_INTERACTION_CK_HPP
