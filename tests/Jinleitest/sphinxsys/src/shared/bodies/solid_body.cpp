#include "solid_body.h"

#include "base_material.h"
#include "base_particles.hpp"
#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
SolidBodyPartForSimbody::
    SolidBodyPartForSimbody(SPHBody &body, Shape &body_part_shape)
    : BodyRegionByParticle(body, body_part_shape),
      rho0_(DynamicCast<Solid>(this, body.getBaseMaterial()).ReferenceDensity()),
      Vol_(base_particles_.getVariableDataByName<Real>("VolumetricMeasure")),
      pos_(base_particles_.getVariableDataByName<Vecd>("Position"))
{
    setMassProperties();
}
//=================================================================================================//
SolidBodyPartForSimbody::SolidBodyPartForSimbody(SPHBody &body, SharedPtr<Shape> shape_ptr)
    : SolidBodyPartForSimbody(body, *shape_ptr.get()) {}
//=================================================================================================//
} // namespace SPH
