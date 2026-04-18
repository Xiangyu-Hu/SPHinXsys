#include "body_part_for_simbody.h"

#include "base_particles.hpp"
namespace SPH
{
//=================================================================================================//
void SolidBodyPartForSimbody::setMassProperties()
{
    Real body_part_volume(0);
    initial_mass_center_ = Vec3d::Zero();
    for (size_t i = 0; i < body_part_particles_.size(); ++i)
    {
        size_t index_i = body_part_particles_[i];
        Vecd particle_position = pos_[index_i];
        Real particle_volume = Vol_[index_i];

        initial_mass_center_ += particle_volume * particle_position;
        body_part_volume += particle_volume;
    }

    initial_mass_center_ /= body_part_volume;

    // computing unit inertia
    Vec3d inertia_moments = Vec3d::Zero();
    Vec3d inertia_products = Vec3d::Zero();
    for (size_t i = 0; i < body_part_particles_.size(); ++i)
    {
        size_t index_i = body_part_particles_[i];
        Vecd particle_position = pos_[index_i];
        Real particle_volume = Vol_[index_i];

        Vec3d displacement = particle_position - initial_mass_center_;

        inertia_moments[0] += particle_volume * (displacement[1] * displacement[1] + displacement[2] * displacement[2]);
        inertia_moments[1] += particle_volume * (displacement[0] * displacement[0] + displacement[2] * displacement[2]);
        inertia_moments[2] += particle_volume * (displacement[0] * displacement[0] + displacement[1] * displacement[1]);
        inertia_products[0] -= particle_volume * displacement[0] * displacement[1];
        inertia_products[1] -= particle_volume * displacement[0] * displacement[2];
        inertia_products[2] -= particle_volume * displacement[1] * displacement[2];
    }
    inertia_moments /= body_part_volume;
    inertia_products /= body_part_volume;

    body_part_mass_properties_ = mass_properties_keeper_.createPtr<SimTK::MassProperties>(
        body_part_volume * rho0_, SimTKVec3(0),
        SimTK::UnitInertia(SimTKVec3(inertia_moments[0], inertia_moments[1], inertia_moments[2]),
                           SimTKVec3(inertia_products[0], inertia_products[1], inertia_products[2])));
}
//=================================================================================================//
} // namespace SPH
