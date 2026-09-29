#include "body_part_for_simbody.h"

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
      rho0_(DynamicCast<Solid>(this, body.getMatterMaterial()).ReferenceDensity()),
      Vol_(base_particles_.getVariableDataByName<Real>("VolumetricMeasure")),
      pos_(base_particles_.getVariableDataByName<Vecd>("Position"))
{
    initialize();
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
    return EigenToSimTK(upgradeToVec3d(initial_mass_center_));
}
//=================================================================================================//
} // namespace SPH
