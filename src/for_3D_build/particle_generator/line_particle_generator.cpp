#include "line_particle_generator.h"

#include "base_body.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
ParticleGenerator<Line>::ParticleGenerator(SPHBody &sph_body)
    : ParticleGenerator<>(sph_body) {}
//=================================================================================================//
void ParticleGenerator<Line>::addLineProperties(
    const Vecd &line_normal, const Vecd &line_binormal, Real thickness, Real width)
{
    line_normal_.push_back(line_normal);
    line_thickness_.push_back(thickness);
    line_binormal_.push_back(line_binormal);
    line_width_.push_back(width);
}
//=================================================================================================//
void ParticleGenerator<Line>::initializeParticleVariables()
{
    ParticleGenerator<>::initializeParticleVariables();
    base_particles_.registerSharedVariableFrom<Vecd>("NormalDirection", line_normal_);
    base_particles_.registerSharedVariableFrom<Vecd>("BinormalDirection", line_binormal_);
    base_particles_.registerSharedVariableFrom<Real>("Thickness", line_thickness_);
    base_particles_.registerSharedVariableFrom<Real>("Width", line_width_);
    base_particles_.addVariableToReload<Vecd>("NormalDirection");
    base_particles_.addVariableToReload<Real>("Thickness");
    base_particles_.addVariableToReload<Vecd>("BinormalDirection");
    base_particles_.addVariableToReload<Real>("Width");
}
//=================================================================================================//
void ParticleGenerator<Line>::initializeParticleVariablesFromReload()
{
    ParticleGenerator<>::initializeParticleVariablesFromReload();
    base_particles_.registerSharedVariableFromReload<Vecd>("NormalDirection");
    base_particles_.registerSharedVariableFromReload<Vecd>("BinormalDirection");
    base_particles_.registerSharedVariableFromReload<Real>("Thickness");
    base_particles_.registerSharedVariableFromReload<Real>("Width");
}
//=================================================================================================//
} // namespace SPH
