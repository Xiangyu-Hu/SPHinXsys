#include "solid_body.h"
#include "base_material.h"
#include "solid_particles.h"
#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
SolidBodyPartForSimbody::
    SolidBodyPartForSimbody(SPHBody &body, Shape &body_part_shape)
    : BodyRegionByParticle(body, body_part_shape),
      solid_body_density_(DynamicCast<Solid>(this, body.base_material_)->ReferenceDensity()),
      solid_particles_(DynamicCast<SolidParticles>(this, &body.getBaseParticles())),
      Vol_(solid_particles_->VolumetricMeasures()),
      pos0_(solid_particles_->pos0_)
{
    setMassProperties();
}
//=================================================================================================//
SolidBodyPartForSimbody::SolidBodyPartForSimbody(SPHBody &body, SharedPtr<Shape> shape_ptr)
    : SolidBodyPartForSimbody(body, *shape_ptr.get()) {}
//=================================================================================================//
} // namespace SPH
