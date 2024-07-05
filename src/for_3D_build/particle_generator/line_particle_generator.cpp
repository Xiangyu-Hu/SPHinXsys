#include "line_particle_generator.h"

#include "base_body.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
ParticleGenerator<LinearParticles>::
    ParticleGenerator(SPHBody &sph_body, LinearParticles &linear_particles)
    : ParticleGenerator<BaseParticles>(sph_body, linear_particles),
      linear_particles_(linear_particles) {}
//=================================================================================================//
void ParticleGenerator<LinearParticles>::addLineProperties(
    const Vecd &line_normal, const Vecd &line_binormal, Real thickness, Real width)
{
    line_normal_.push_back(line_normal);
    line_thickness_.push_back(thickness);
    line_binormal_.push_back(line_binormal);
    line_width_.push_back(width);
}
//=================================================================================================//
void ParticleGenerator<LinearParticles>::initializeParticleVariables()
{
    ParticleGenerator<BaseParticles>::initializeParticleVariables();
    linear_particles_.registerSurfaceProperties(line_normal_, line_thickness_);
    linear_particles_.registerLineProperties(line_binormal_, line_width_);
}
//=================================================================================================//
void ParticleGenerator<LinearParticles>::initializeParticleVariablesFromReload()
{
    ParticleGenerator<BaseParticles>::initializeParticleVariablesFromReload();
    linear_particles_.registerSurfacePropertiesFromReload();
    linear_particles_.registerLinePropertiesFromReload();
}
//=================================================================================================//
} // namespace SPH
