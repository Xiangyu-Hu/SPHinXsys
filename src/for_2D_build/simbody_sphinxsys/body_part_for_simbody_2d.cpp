#include "body_part_for_simbody.h"

#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
void SolidBodyPartForSimbody::initialize()
{
    Real body_part_volume(0);
    Vecd mass_center = Vecd::Zero();
    for (size_t i = 0; i < body_part_particles_.size(); ++i)
    {
        size_t index_i = body_part_particles_[i];
        Real particle_volume = Vol_[index_i];
        mass_center += particle_volume * pos_[index_i];
        body_part_volume += particle_volume;
    }

    mass_center /= body_part_volume;
    initial_mass_center_ = mass_center;

    // computing unit inertia
    Real Ixx = 0.0, Iyy = 0.0, Ixy = 0.0; // Ixz = Iyz = 0 in 2D
    for (size_t i = 0; i < body_part_particles_.size(); ++i)
    {
        size_t index_i = body_part_particles_[i];
        Vecd r = pos_[index_i] - mass_center; // relative position vector
        Real x = r[0];
        Real y = r[1];
        // z = 0 for planar case
        Real particle_volume = Vol_[index_i];

        Ixx += particle_volume * y * y;
        Iyy += particle_volume * x * x;
        Ixy -= particle_volume * x * y; // note the negative sign
    }
    Ixx /= body_part_volume;
    Iyy /= body_part_volume;
    Ixy /= body_part_volume;

    Real Izz = Ixx + Iyy; // For 2D, Izz = Ixx + Iyy
    SimTK::UnitInertia unit_inertia(
        Ixx, Iyy, Izz, // diagonal moments
        Ixy, 0.0, 0.0  // products of inertia
    );

    // Zero center of mass due to local frame
    body_part_mass_properties_ = mass_properties_keeper_.createPtr<SimTK::MassProperties>(
        body_part_volume * rho0_, SimTKVec3(0.0), unit_inertia);
}
//=================================================================================================//
} // namespace SPH