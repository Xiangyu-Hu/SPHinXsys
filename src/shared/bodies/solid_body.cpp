#include "solid_body.h"
#include "base_material.h"
#include "solid_particles.h"
#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
SolidBodyPartForSimbody::
    SolidBodyPartForSimbody(SPHBody &body, SharedPtr<Shape> shape_ptr)
    : BodyRegionByParticle(body, shape_ptr),
      solid_particles_(DynamicCast<SolidParticles>(this, &body.getBaseParticles())),
      solid_body_density_(solid_particles_->getBaseMaterial().ReferenceDensity())
{
    setMassProperties();
}
//=================================================================================================//
} // namespace SPH
