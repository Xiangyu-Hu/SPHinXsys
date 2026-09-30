#include "body_part_for_simbody.h"

#include "adaptation.h"
#include "base_body.h"
#include "base_material.h"
#include "base_particles.hpp"
#include "vector_functions.h"

namespace SPH
{
//=================================================================================================//
SolidBodyPartForSimbody::
    SolidBodyPartForSimbody(SPHBody &body, Shape &body_part_shape)
    : BodyRegionByParticle(body, body_part_shape),
      min_spacing_sqr_(pow(body.getSPHAdaptation().MinimumSpacing(), 2)),
      rho0_(DynamicCast<Solid>(this, body.getMatterMaterial()).ReferenceDensity()),
      Vol_(base_particles_.getVariableDataByName<Real>("VolumetricMeasure")),
      pos_(base_particles_.getVariableDataByName<Vecd>("Position"))
{
    Real body_part_volume(0);
    Vec3d mass_center = Vec3d::Zero();
    for (size_t i = 0; i < body_part_particles_.size(); ++i)
    {
        size_t index_i = body_part_particles_[i];
        Vec3d particle_position = upgradeToVec3d(pos_[index_i]);
        Real particle_volume = Vol_[index_i];

        mass_center += particle_volume * particle_position;
        body_part_volume += particle_volume;
    }
    mass_center /= body_part_volume;

    // computing unit inertia
    Vec3d inertia_moments = Vec3d::Zero();
    Vec3d inertia_products = Vec3d::Zero();
    for (size_t i = 0; i < body_part_particles_.size(); ++i)
    {
        size_t index_i = body_part_particles_[i];
        Vec3d particle_position = upgradeToVec3d(pos_[index_i]);
        Real particle_volume = Vol_[index_i];

        Vec3d displacement = particle_position - mass_center;

        inertia_moments[0] += particle_volume * (displacement[1] * displacement[1] + displacement[2] * displacement[2]);
        inertia_moments[1] += particle_volume * (displacement[0] * displacement[0] + displacement[2] * displacement[2]);
        inertia_moments[2] += particle_volume * (displacement[0] * displacement[0] + displacement[1] * displacement[1]);
        inertia_products[0] -= particle_volume * displacement[0] * displacement[1];
        inertia_products[1] -= particle_volume * displacement[0] * displacement[2];
        inertia_products[2] -= particle_volume * displacement[1] * displacement[2];
    }
    inertia_moments /= body_part_volume;
    inertia_products /= body_part_volume;

    initial_mass_center_ = EigenToSimTK(mass_center);
    if constexpr (Dimensions == 2)
    {
        inertia_moments[2] -= min_spacing_sqr_ / 6.0; // assume a small thickness
    }
    SimTK::UnitInertia unit_inertia(
        EigenToSimTK(inertia_moments), EigenToSimTK(inertia_products));
    // Zero center of mass due to local frame
    body_part_mass_properties_ = mass_properties_keeper_.createPtr<SimTK::MassProperties>(
        body_part_volume * rho0_, SimTK::Vec3(Real(0)), unit_inertia);
}
//=================================================================================================//
SolidBodyPartForSimbody::SolidBodyPartForSimbody(SPHBody &body, SharedPtr<Shape> shape_ptr)
    : SolidBodyPartForSimbody(body, *shape_ptr.get()) {}
//=================================================================================================//
SimTK::MassProperties &SolidBodyPartForSimbody::getSimTKMassProperties() const
{
    return *body_part_mass_properties_;
}
//=================================================================================================//
SimTK::Vec3 SolidBodyPartForSimbody::getSimTKMassCenter() const
{
    return initial_mass_center_;
}
//=================================================================================================//
} // namespace SPH